# hash:sha256:ff6a38592f47ddceee71bef0334dbb7c875eac796288b887dfbb34d3d81631f5
FROM registry.codeocean.com/codeocean/matlab:2021b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

COPY postInstall /
RUN /postInstall
