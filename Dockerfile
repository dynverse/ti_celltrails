FROM dynverse/dynwrapr:v0.1.0

RUN R -e 'devtools::install_github("dcellwanger/CellTrails")'

COPY definition.yml run.R example.sh /code/

ENTRYPOINT ["/code/run.R"]
