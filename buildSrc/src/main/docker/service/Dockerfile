
# NOTE: the build context for this Dockerfile is the project root

# use JDK 17 as the base image, since that's probably the heaviest dependency here
# image based on Ubuntu 20.04 LTS (Focal Fossa)
# https://hub.docker.com/_/eclipse-temurin/
FROM eclipse-temurin:17.0.2_8-jdk-focal

ARG VERSIONS

# install s6-overlay as the init system
# https://github.com/just-containers/s6-overlay
ADD https://github.com/just-containers/s6-overlay/releases/download/v2.2.0.1/s6-overlay-amd64-installer /tmp/
RUN chmod +x /tmp/s6-overlay-amd64-installer && /tmp/s6-overlay-amd64-installer /
ENTRYPOINT ["/init"]

# set up directories
ARG servicedir=buildSrc/src/main/docker/service
WORKDIR /opt/osprey

# add the service versions
ADD build/docker/versions versions
LABEL "osprey.service.versions"="$VERSIONS"

# install libraries needed by AmberTools
RUN apt-get update \
	&& apt-get install -y libgfortran4

# install s6 services for all the API versions
COPY $servicedir/install-osprey.sh ./
COPY $servicedir/finish.sh ./
RUN chmod +x install-osprey.sh \
	&& chmod +x finish.sh \
	&& ./install-osprey.sh

# install Caddy 2
# https://caddyserver.com/v2
# https://github.com/caddyserver/caddy-docker/blob/master/Dockerfile.tmpl
ENV CADDY_VERSION 2.4.3
ADD https://github.com/caddyserver/caddy/releases/download/v$CADDY_VERSION/caddy_${CADDY_VERSION}_linux_amd64.tar.gz /tmp/caddy.tar.gz
RUN tar x -z -f /tmp/caddy.tar.gz -C /usr/bin caddy \
	&& rm -f /tmp/caddy.tar.gz \
	&& chmod +x /usr/bin/caddy \
	&& caddy version

# configure caddy
COPY $servicedir/Caddyfile ./

# install s6 service for caddy
# if caddy ever dies, tell s6 to exit and not to restart it
# ideally, the whole container should die if any one service dies
COPY $servicedir/start-caddy.sh /etc/services.d/caddy/run
RUN chmod +x /etc/services.d/caddy/run \
	&& cp finish.sh /etc/services.d/caddy/finish

# create a place the service can store data
RUN mkdir -p /home/osprey \
	&& chown www-data: /home/osprey \
	&& ls -al /home/osprey
VOLUME /home/osprey
