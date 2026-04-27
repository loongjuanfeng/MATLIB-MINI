# MATLIB-MINI

C++23 header-only linear algebra library.

## Build

```bash
# Release
cmake --preset default && ninja -C build

# Run demo
./build/main
```

## Test

```bash
make test         # Release
make test-debug   # Debug
```

Requires C++23 and TBB (fetched automatically if not installed).
