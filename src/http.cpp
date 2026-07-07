#include "http.hpp"

#include <curl/curl.h>

#include <chrono>
#include <mutex>
#include <stdexcept>
#include <thread>

namespace sed::http {

namespace {

// One-time curl_global_init, torn down at process exit.
struct GlobalInit {
  GlobalInit() { curl_global_init(CURL_GLOBAL_DEFAULT); }
  ~GlobalInit() { curl_global_cleanup(); }
};
void ensure_global() { static GlobalInit g; }

std::mutex g_throttle_mutex;
int g_throttle_ms = 0;
std::chrono::steady_clock::time_point g_last_request{};

// Block until at least g_throttle_ms have elapsed since the previous request
// start, then record this request's start time.
void throttle_gate() {
  std::lock_guard<std::mutex> lock(g_throttle_mutex);
  if (g_throttle_ms > 0) {
    auto now = std::chrono::steady_clock::now();
    auto wait = std::chrono::milliseconds(g_throttle_ms) - (now - g_last_request);
    if (wait.count() > 0) std::this_thread::sleep_for(wait);
  }
  g_last_request = std::chrono::steady_clock::now();
}

size_t write_cb(char* ptr, size_t size, size_t nmemb, void* userdata) {
  auto* out = static_cast<std::string*>(userdata);
  out->append(ptr, size * nmemb);
  return size * nmemb;
}

// Shared setup for GET/POST. Returns the performed Response or throws.
// max_tls12: some legacy archive servers (e.g. archive.stsci.edu/pub) fail the
// default TLS 1.3 handshake; on CURLE_SSL_CONNECT_ERROR we retry once capped
// at TLS 1.2 (wget/GnuTLS negotiates these hosts fine, so ISIS never saw it).
Response perform(const std::string& url, const Options& opt, bool is_post,
                 const std::string& body, const std::string& content_type,
                 bool max_tls12 = false) {
  ensure_global();
  throttle_gate();

  CURL* curl = curl_easy_init();
  if (!curl) throw std::runtime_error("curl_easy_init failed");

  Response resp;
  struct curl_slist* header_list = nullptr;

  curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_cb);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, &resp.body);
  curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, (long)opt.connect_timeout_s);
  curl_easy_setopt(curl, CURLOPT_TIMEOUT, (long)opt.timeout_s);
  curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, opt.follow_redirects ? 1L : 0L);
  curl_easy_setopt(curl, CURLOPT_USERAGENT, "sedfit/2.0 (+libcurl)");
  curl_easy_setopt(curl, CURLOPT_ACCEPT_ENCODING, "");  // accept gzip/deflate

  std::string all_headers;
  for (const auto& h : opt.headers)
    header_list = curl_slist_append(header_list, h.c_str());

  if (is_post) {
    std::string ctype = "Content-Type: " + content_type;
    header_list = curl_slist_append(header_list, ctype.c_str());
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE, (long)body.size());
  }
  if (header_list) curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header_list);
  if (max_tls12)
    curl_easy_setopt(curl, CURLOPT_SSLVERSION,
                     CURL_SSLVERSION_TLSv1_0 | CURL_SSLVERSION_MAX_TLSv1_2);

  CURLcode rc = curl_easy_perform(curl);
  if (rc == CURLE_OK)
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &resp.status);

  if (header_list) curl_slist_free_all(header_list);
  curl_easy_cleanup(curl);

  if (rc == CURLE_SSL_CONNECT_ERROR && !max_tls12)
    return perform(url, opt, is_post, body, content_type, /*max_tls12=*/true);
  if (rc != CURLE_OK)
    throw std::runtime_error(std::string("HTTP request failed: ") +
                             curl_easy_strerror(rc) + " (" + url + ")");
  return resp;
}

}  // namespace

void set_throttle_ms(int ms) {
  std::lock_guard<std::mutex> lock(g_throttle_mutex);
  g_throttle_ms = ms;
}

Response get(const std::string& url, const Options& opt) {
  return perform(url, opt, /*is_post=*/false, "", "");
}

Response post(const std::string& url, const std::string& body,
              const std::string& content_type, const Options& opt) {
  return perform(url, opt, /*is_post=*/true, body, content_type);
}

bool probe(const std::string& url, int connect_timeout_s, int timeout_s) {
  try {
    Options opt;
    opt.connect_timeout_s = connect_timeout_s;
    opt.timeout_s = timeout_s;
    Response r = get(url, opt);
    return r.status != 0;  // any HTTP answer means reachable
  } catch (...) {
    return false;
  }
}

std::string url_encode(const std::string& s) {
  ensure_global();
  CURL* curl = curl_easy_init();
  if (!curl) throw std::runtime_error("curl_easy_init failed");
  char* enc = curl_easy_escape(curl, s.c_str(), (int)s.size());
  std::string out = enc ? enc : "";
  if (enc) curl_free(enc);
  curl_easy_cleanup(curl);
  return out;
}

}  // namespace sed::http
