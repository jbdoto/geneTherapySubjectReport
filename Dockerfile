FROM  continuumio/anaconda3

# Copy project files into container
# (not globbing all files to prevent IntelliJ IDE files from breaking Docker caching):
COPY ./report.R /report/report.R

COPY conda_spec.txt /report/conda_spec.txt

WORKDIR /report
RUN apt-get update -y  && \
    apt-get install -y libcurl4-openssl-dev \
    libmariadb-dev \
    gcc \
    g++ \
    make  && \
    apt-get clean

ENV C_INCLUDE_PATH=/usr/lib/gcc/x86_64-linux-gnu/
RUN ln -s /usr/lib/x86_64-linux-gnu/libmariadb.so /usr/lib/x86_64-linux-gnu/libmysqlclient.so.18

RUN conda create --name report --file conda_spec.txt


# Make RUN commands use the new environment:
# https://pythonspeed.com/articles/activate-conda-dockerfile/
SHELL ["conda", "run", "-n", "report", "/bin/bash", "-c"]

# Install non-Conda dependencies:
COPY installPackages.R /report/installPackages.R
RUN Rscript installPackages.R
WORKDIR /scratch
ENTRYPOINT ["conda", "run", "-n", "report", "Rscript", "/report/report.R"]
