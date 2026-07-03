// Thin libcurl wrapper for the query layer: GET/POST with timeouts, a global
// init guard, a connectivity probe, and a global min-interval throttle so bulk
// multi-star runs do not hammer the archives.
#pragma once

#include <string>
#include <vector>

namespace sed::http {

struct Response {
  long status = 0;        // HTTP status code (0 if no response)
  std::string body;       // response body
  bool ok() const { return status >= 200 && status < 300; }
};

// Options for a single request.
struct Options {
  int connect_timeout_s = 15;  // TCP connect timeout
  int timeout_s = 120;         // total transfer timeout
  bool follow_redirects = true;
  std::vector<std::string> headers;  // extra "Key: Value" request headers
};

// Set a global minimum interval (milliseconds) enforced between the *start* of
// consecutive requests across all threads. 0 disables throttling. Default 0.
void set_throttle_ms(int ms);

// HTTP GET. Throws std::runtime_error on transport failure (DNS, connect,
// timeout); a non-2xx HTTP status is returned in Response, not thrown.
Response get(const std::string& url, const Options& opt = {});

// HTTP POST with a raw body and an explicit content type.
Response post(const std::string& url, const std::string& body,
              const std::string& content_type = "application/x-www-form-urlencoded",
              const Options& opt = {});

// Reachability probe: returns true if the URL answers within the timeouts
// (mirrors the S-Lang `curl -fs --connect-timeout N -m M ... ; echo $?`==0
// gate). Never throws.
bool probe(const std::string& url, int connect_timeout_s = 4, int timeout_s = 6);

// Percent-encode a string for use in a URL query component.
std::string url_encode(const std::string& s);

}  // namespace sed::http
