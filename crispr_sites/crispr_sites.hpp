// Look for 20-mers at PAM sites.  Including NGG or CCN, k=23.
constexpr auto k = 23;

// input will be read STRIDE_SIZE at a time
constexpr auto STRIDE_SIZE = 32 * 1024 * 1024;

// to scan for k-mers, consecutive read windows must overlap by k-1 characters
constexpr auto BUFFER_SIZE = STRIDE_SIZE + k - 1;
