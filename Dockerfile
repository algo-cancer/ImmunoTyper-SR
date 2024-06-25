FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /app

# Install git
RUN apt-get update && apt-get install -y git && apt-get clean

# Clone the repository
RUN git clone https://github.com/michael-ford/ImmunoTyper-SR.git

# Create the conda environment and install required packages
RUN source /opt/conda/etc/profile.d/conda.sh
RUN conda config --prepend channels conda-forge
RUN conda config --append channels bioconda
RUN conda install -y -c bioconda python=3.8 bwa bowtie \
    && conda install -y -c gurobi gurobi && conda install -y -c conda-forge samtools

# Install the ImmunoTyper-SR package
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && pip install /app/ImmunoTyper-SR"

# Set the Gurobi license environment variable
ENV GRB_LICENSE_FILE=/opt/gurobi/gurobi.lic

# Copy the entrypoint script into the container
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["--help"]