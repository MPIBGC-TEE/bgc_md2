FROM python:3.9.16
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
WORKDIR ${HOME}
USER ${NB_USER}
ENV PATH="${HOME}/.local/bin:${PATH}"
RUN python -m pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r req.txt notebook jupyterlab
RUN pwd;ls -la
