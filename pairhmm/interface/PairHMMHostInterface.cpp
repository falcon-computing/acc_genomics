#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <vector>

#include "ksight/tools.h"

#include "PairHMMHostInterface.h"

// Encode a scalar value to serialized data
template <typename T>
static inline void putT(void* buf, uint64_t & buf_idx, T value) {
  int len = sizeof(value);
  memcpy((char*)buf + buf_idx, reinterpret_cast<char*>(&value), len);
  buf_idx += len;
}

// Decode a scalar value from serialized data
template <typename T>
static inline void getT(const void* buf, uint64_t & buf_idx, T &value) {
  int len = sizeof(value);
  memcpy(reinterpret_cast<char*>(&value), (char*)buf + buf_idx, len);
  buf_idx += len;
}

// Store a string with its length to serialized data
static inline void putStr(void* buf, uint64_t & buf_idx, const char* str, size_t len) {
  memcpy((char*)buf + buf_idx, str, len);
  buf_idx += len;
}

// Retrieve a string from serialized data
static inline void getStr(const void* buf, uint64_t & buf_idx, char* &dst, size_t len) {
  if (len > 0) {
    dst = (char*)malloc(len+1);
    memcpy(dst, (char*)buf + buf_idx,  len);
    dst[len] = '\0'; // add trailing \0 for compatibility
    buf_idx += len;
  }
}

// Store a string with its length to serialized data
static inline void putStr(void* buf, uint64_t & buf_idx, const std::string &str) {
  memcpy((char*)buf + buf_idx, str.c_str(), str.length());
  buf_idx += str.length();
}

// Retrieve a string from serialized data
static inline void getStr(const void* buf, uint64_t & buf_idx, std::string &data, uint64_t length) {
  if (length > 0) {
    data = std::string((char*)buf + buf_idx, length);
    buf_idx += length;
  }
}

// Encode a scalar value to serialized data
template <typename T>
static inline void putT(std::stringstream &ss, T value) {
  ss.write(reinterpret_cast<char*>(&value), sizeof(value));
}

// Decode a scalar value from serialized data
template <typename T>
static inline void getT(std::stringstream &ss, T &value) {
  ss.read(reinterpret_cast<char*>(&value), sizeof(value));
}

// Store a string with its length to serialized data
static inline void putStr(std::stringstream &ss, const std::string &str) {
  size_t length = str.size();
  //putT(ss, length);
  ss.write(str.c_str(), length);
  //printf("put %d\n", length);
}

// Retrieve a string from serialized data
static inline void getStr(std::stringstream &ss, 
    uint64_t length, 
    std::string & data) {
  data.resize(length);
  ss.read(&data[0], length);
}

// Store a string with its length to serialized data
static inline void putStr(std::stringstream &ss, const char* str, size_t len) {
  ss.write(str, len);
}

// Retrieve a string from serialized data
static inline void getStr(std::stringstream &ss, char* &dst, size_t len) {
  if (len > 0) {
    dst = (char*)malloc(len+1);
    ss.read(dst, len);
    dst[len] = '\0'; // add trailing \0 for compatibility
  }
}

#if 0
std::string serialize(const std::vector<Read> & reads)
{
  PLACE_TIMER;
  std::stringstream ss;
  int num = reads.size();
  putT(ss, num);
  for (int i = 0; i < reads.size(); i++) {
    int len = reads[i].bases.size(); 

    putT(ss, len);
    putStr(ss, reads[i].bases);
    putStr(ss, reads[i]._q);
    putStr(ss, reads[i]._i);
    putStr(ss, reads[i]._d);
    putStr(ss, reads[i]._c);
  }
  return ss.str();
}

std::string serialize(const std::vector<Hap> & haps)
{
  PLACE_TIMER;
  std::stringstream ss;
  int num = haps.size();
  putT(ss, num);
  for (int i = 0; i < haps.size(); i++) {
    int len = haps[i].bases.size(); 

    putT(ss, len);
    putStr(ss, haps[i].bases);
  }
  return ss.str();
}

void deserialize(const std::string & data,
    std::vector<Read> & reads) {
  PLACE_TIMER;
  
  std::stringstream ss(data);

  int num = 0;
  getT(ss, num);
  reads.resize(num);

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(ss, len);

    getStr(ss, len, reads[i].bases);
    getStr(ss, len, reads[i]._q);
    getStr(ss, len, reads[i]._i);
    getStr(ss, len, reads[i]._d);
    getStr(ss, len, reads[i]._c);
  }
}

void deserialize(const std::string & data,
    std::vector<Hap> & haps) {
  PLACE_TIMER;
  
  std::stringstream ss(data);

  int num = 0;
  getT(ss, num);
  haps.resize(num);

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(ss, len);
    
    getStr(ss, len, haps[i].bases);
  }
}
#endif

uint64_t serialize(void* buf, const read_t * reads, int num) {
  PLACE_TIMER1("reads");

  uint64_t buf_idx = 0;
  putT(buf, buf_idx, num);

  for (int i = 0; i < num; i++) {
    int len = reads[i].len; 
    putT(buf, buf_idx, len);

    putStr(buf, buf_idx, reads[i]._b, len);
    putStr(buf, buf_idx, reads[i]._q, len);
    putStr(buf, buf_idx, reads[i]._i, len);
    putStr(buf, buf_idx, reads[i]._d, len);
    putStr(buf, buf_idx, reads[i]._c, len);
  }
  return buf_idx;
}

uint64_t serialize(void* buf, const hap_t * haps, int num) {
  PLACE_TIMER1("haps");

  uint64_t buf_idx = 0;
  putT(buf, buf_idx, num);

  for (int i = 0; i < num; i++) {
    int len = haps[i].len;
    putT(buf, buf_idx, len);
    putStr(buf, buf_idx, haps[i]._b, len);
  }
  return buf_idx;
}

int deserialize(
    const void* buf, 
    read_t* &reads) 
{
  PLACE_TIMER1("reads");

  uint64_t buf_idx = 0;
  int num = 0;
  getT(buf, buf_idx, num);

  reads = (read_t*)malloc(num*sizeof(read_t));

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(buf, buf_idx, len);
    reads[i].len = len;

    getStr(buf, buf_idx, reads[i]._b, len);
    getStr(buf, buf_idx, reads[i]._q, len);
    getStr(buf, buf_idx, reads[i]._i, len);
    getStr(buf, buf_idx, reads[i]._d, len);
    getStr(buf, buf_idx, reads[i]._c, len);
  }
  return num;
}

int deserialize(
    const void* buf, 
    hap_t* &haps) 
{
  PLACE_TIMER1("haps");
  
  uint64_t buf_idx = 0;
  int num = 0;
  getT(buf, buf_idx, num);

  haps = (hap_t*)malloc(num*sizeof(hap_t));

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(buf, buf_idx, len);
    haps[i].len = len;

    getStr(buf, buf_idx, haps[i]._b, len);
  }
  
  return num;
}

std::string serialize(const read_t * reads, int num) {
  PLACE_TIMER1("reads");

  std::stringstream ss;
  putT(ss, num);

  for (int i = 0; i < num; i++) {
    int len = reads[i].len; 
    putT(ss, len);

    putStr(ss, reads[i]._b, len);
    putStr(ss, reads[i]._q, len);
    putStr(ss, reads[i]._i, len);
    putStr(ss, reads[i]._d, len);
    putStr(ss, reads[i]._c, len);
  }

  return ss.str();
}

std::string serialize(const hap_t * haps, int num) {
  PLACE_TIMER1("haps");

  std::stringstream ss;
  putT(ss, num);

  for (int i = 0; i < num; i++) {
    int len = haps[i].len;
    putT(ss, len);

    putStr(ss, haps[i]._b, len);
  }

  return ss.str();
}

int deserialize(const std::string & data,
    read_t* &reads) 
{
  PLACE_TIMER1("reads");
  
  std::stringstream ss(data);

  int num = 0;
  getT(ss, num);

  reads = (read_t*)malloc(num*sizeof(read_t));

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(ss, len);
    reads[i].len = len;

    getStr(ss, reads[i]._b, len);
    getStr(ss, reads[i]._q, len);
    getStr(ss, reads[i]._i, len);
    getStr(ss, reads[i]._d, len);
    getStr(ss, reads[i]._c, len);
  }
  return num;
}

int deserialize(const std::string & data,
    hap_t* &haps) {
  PLACE_TIMER1("haps");
  
  std::stringstream ss(data);

  int num = 0;
  getT(ss, num);

  haps = (hap_t*)malloc(num*sizeof(hap_t));

  for (int i = 0; i < num; i++) {
    int len = 0;
    getT(ss, len);
    haps[i].len = len;
    getStr(ss, haps[i]._b, len);
  }
  
  return num;
}
