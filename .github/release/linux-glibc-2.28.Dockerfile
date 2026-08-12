FROM rockylinux/rockylinux:8.10

ARG NIM_VERSION=2.2.10
ARG NIM_LINUX_X64_SHA256=0a3a38752e97e9d44aa479b3a7b37336dfe0176daf22ee5b5218ad0991ecd211

RUN dnf -y install findutils \
      binutils \
      ca-certificates \
      curl \
      file \
      findutils \
      gcc \
      glibc-devel \
      make \
      python39 \
      tar \
      xz \
    && dnf clean all \
    && rm -rf /var/cache/dnf

RUN set -eux; \
    test "$(getconf GNU_LIBC_VERSION | awk '{print $2}')" = "2.28"; \
    archive="nim-${NIM_VERSION}-linux_x64.tar.xz"; \
    curl --fail --location --show-error --silent \
      --retry 5 --retry-delay 2 --retry-connrefused \
      --output "/tmp/${archive}" \
      "https://nim-lang.org/download/${archive}"; \
    printf '%s  %s\n' "${NIM_LINUX_X64_SHA256}" "/tmp/${archive}" | sha256sum -c -; \
    tar --extract --xz --no-same-owner --file "/tmp/${archive}" --directory /opt; \
    test -x "/opt/nim-${NIM_VERSION}/bin/nim"; \
    rm -f "/tmp/${archive}"

ENV PATH="/opt/nim-2.2.10/bin:${PATH}"
WORKDIR /work
