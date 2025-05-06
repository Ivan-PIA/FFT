#!/bin/bash

set -e 

MODE=${1:-run} 

BUILD_DIR=build

if [ "$MODE" = "test" ]; then
    echo "--- Building with tests enabled..."
    cmake -B $BUILD_DIR -DBUILD_TESTS=ON
    cmake --build $BUILD_DIR
    $BUILD_DIR/test_fft
else
    echo "--- Building main program..."
    cmake -B $BUILD_DIR
    cmake --build $BUILD_DIR
    $BUILD_DIR/fft_main
fi