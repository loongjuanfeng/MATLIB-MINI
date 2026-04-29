.PHONY: all configure build run test test-debug clean

BUILD_DIR   := build
BUILD_DEBUG := build-debug

# --- Build ---
all: build

configure:
	cmake --preset default

build: configure
	ninja -C $(BUILD_DIR)

run: build
	$(BUILD_DIR)/main

# --- Test ---
test: configure
	ninja -C $(BUILD_DIR)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

test-debug:
	cmake --preset debug
	ninja -C $(BUILD_DEBUG)
	ctest --test-dir $(BUILD_DEBUG) --output-on-failure

# --- Clean ---
clean:
	rm -rf $(BUILD_DIR) $(BUILD_DEBUG) .cache
