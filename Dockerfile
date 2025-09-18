FROM ubuntu:22.04

ARG STAR_VERSION=2.7.11b
ENV STAR_VERSION=${STAR_VERSION}
ARG BOWTIE2_VERSION=2.5.4
ENV BOWTIE2_VERSION=${BOWTIE2_VERSION}

# Avoid interactive prompts and set a default timezone
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    build-essential \
    g++ \
    make \
    perl \
    python3 \
    r-base \
    tzdata \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    bowtie \
    xxd \
    samtools \
    tabix \
    unzip \
    wget \
    && rm -rf /var/lib/apt/lists/* \
    && ln -fs /usr/share/zoneinfo/$TZ /etc/localtime \
    && dpkg-reconfigure -f noninteractive tzdata

# Build STAR and Bowtie2 with AVX optimizations
RUN set -eux; \
    if grep -qi 'avx512' /proc/cpuinfo; then \
        SIMD_FLAGS="-mavx512f -mavx512bw -mavx512dq -mavx512cd -mavx512vl -mavx2"; \
    elif grep -qi 'avx2' /proc/cpuinfo; then \
        SIMD_FLAGS="-mavx2"; \
    else \
        echo "CPU does not support AVX2 or AVX-512 instructions" >&2; \
        exit 1; \
    fi; \
    echo "Using SIMD flags: ${SIMD_FLAGS}"; \
    TMPDIR="$(mktemp -d)"; \
    cd "${TMPDIR}"; \
    wget -qO star.tar.gz "https://github.com/alexdobin/STAR/archive/refs/tags/${STAR_VERSION}.tar.gz"; \
    tar -xzf star.tar.gz; \
    cd STAR-${STAR_VERSION}/source; \
    make -j"$(nproc)" STAR CXXFLAGSextra="${SIMD_FLAGS}"; \
    make -j"$(nproc)" STARlong CXXFLAGSextra="${SIMD_FLAGS}"; \
    install -m 0755 STAR /usr/local/bin/STAR; \
    install -m 0755 STARlong /usr/local/bin/STARlong; \
    cd "${TMPDIR}"; \
    wget -qO bowtie2.tar.gz "https://github.com/BenLangmead/bowtie2/archive/refs/tags/v${BOWTIE2_VERSION}.tar.gz"; \
    tar -xzf bowtie2.tar.gz; \
    cd bowtie2-${BOWTIE2_VERSION}; \
    make -j"$(nproc)" SSE_FLAG="${SIMD_FLAGS} -faligned-new -DSSE_AVX2"; \
    make install PREFIX=/usr/local SSE_FLAG="${SIMD_FLAGS} -faligned-new -DSSE_AVX2"; \
    cd /; \
    rm -rf "${TMPDIR}"

# Set working directory
WORKDIR /opt/rsem

# Copy source
COPY . /opt/rsem

# Build and install RSEM
RUN make && make install

ENV PATH="/usr/local/bin:/usr/bin:${PATH}"

CMD ["bash"]
