FROM --platform=linux/amd64 ubuntu:22.04

# mpboot requires x86 SSE3 instructions.
# On Apple Silicon (M1/M2/M3), enable Rosetta in Docker Desktop:
#   Settings → General → "Use Rosetta for x86_64/amd64 emulation on Apple Silicon"

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy source and build
COPY . /app/source/
RUN mkdir -p /app/build && cd /app/build \
    && cmake /app/source -DCMAKE_BUILD_TYPE=Release -DIQTREE_FLAGS="sse" \
    && make -j$(nproc)

RUN chmod +x /app/source/scripts/run_spr_benchmark.sh

ENV PATH="/app/build:${PATH}"
ENV MPBOOT_BUILD_DIR="/app/build"
ENV MPBOOT_DATA_DIR="/app/placement_genbank_data"

ENTRYPOINT ["/app/source/scripts/run_spr_benchmark.sh"]
CMD []
