include platform.make

all: test

test: bin/shade
	@bin/shade

bin/shade: src/shade.cpp build/fragment.glsl.h
	@echo "C++ $<"
	@mkdir -p bin
	@$(cxx) -std=c++14 -Os -I. $< -o $@ $(ldflags)

build/fragment.glsl.h: fragment.glsl scripts/headerify.sh
	@echo "HDR $<"
	@mkdir -p build
	@scripts/headerify.sh $< > $@

clean:
	@echo "CLEAN"
	@rm -fr build

.PHONY: all test clean
