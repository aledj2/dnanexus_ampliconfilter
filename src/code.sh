#!/bin/bash


# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

mark-section "download inputs"

# Download inputs
dx-download-all-inputs --parallel

# Create output/input folders
mkdir -p out/discarded_BAM/output/ out/discarded_BAM_BAI/output/ out/clipped_BAM_BAI/output/ out/clipped_BAM/output/ out/metrics/QC/ input_files


# Store the DNA Nexus API key. Grants the script access to DNAnexus resources    
API_KEY=$(dx cat project-FQqXfYQ0Z0gqx7XG9Z2b4K43:mokaguys_nexus_auth_key)

# download the reference genome from 001_Tools...
dx download project-ByfFPz00jy1fk6PjpZ95F27J:file-ByYgX700b80gf4ZY1GxvF3Jv --auth ${API_KEY}
# extract reference genome
tar -xf hs37d5.fasta-index.tar.gz


# SET VARIABLES
#
# Set up options - adding each optional argument to variable 
#
opts=" --clipping ${clipping} "
# add binary flags
if [ "$super" == "true" ]; then
  opts="$opts --super"
fi
if [ "$mask" == "true" ]; then
  opts="$opts --mask"
fi
# add arguments with values
if [ "$primerdistance" != "" ]; then
  opts="$opts --primerdistance ${primerdistance}"
fi
if [ "$tolerance" != "" ]; then
  opts="$opts --tolerance ${tolerance}"
fi
if [ "$maxbuffer" != "" ]; then
  opts="$opts --maxbuffer ${maxbuffer}"
fi
if [ "$extratrim" != "" ]; then
  opts="$opts --extratrim ${extratrim}"
fi

# capture the samplename, exlcuding anything after the first full stop 
# NB the expected input is the refined.bam so will have to add in .refined when manually renaming outputs below
samplename=${BAM_name%%.*}

# move the inputs into the dir which will be mounted in the docker image.
mv genome.fa input_files/
mv $BEDPE_path input_files/
mv $BAM_path input_files/

# download the docker file from 001
dx download project-ByfFPz00jy1fk6PjpZ95F27J:file-FqZvbxj0jy1vvFxvJP6pqKq7 --auth ${API_KEY}
# decompress docker image
tar -xf AmpliconFilter_v1.0.tar.gz
# load docker image
docker load  --input AmpliconFilter_v1.0.tar
# docker run. mount input directory as /sandbox - all inputs were moved there earlier, and all outputs will be saved there
# pass the BEDPE, BAM and reference genome files as inputs, along with $opts string
# name output BAMS with generic names as these will be renamed when sorting and indexing.
# output metrics will be named using the samplename
docker run -v /home/dnanexus/input_files:/sandbox ampliconfilter:1.0 /sandbox/$BEDPE_name -i /sandbox/$BAM_name -g /sandbox/genome.fa $opts -o /sandbox/primerclipped.bam  -d /sandbox/discarded.bam  -m /sandbox/$samplename.refined.primerclipped.metrics

# Downstream tools *may* need indexed BAMs. To index BAMs first need to sort.
# sort BAM and give prefix sorted (created sorted.bam)
samtools sort input_files/primerclipped.bam sorted
# move to output folder and and rename 
mv sorted.bam out/clipped_BAM/output/$samplename.refined.primerclipped.bam 
# index and then move index into seperate output folder
samtools index out/clipped_BAM/output/$samplename.refined.primerclipped.bam 
mv out/clipped_BAM/output/$samplename.refined.primerclipped.bam.bai out/clipped_BAM_BAI/output/$samplename.refined.primerclipped.bam.bai

#repeat for discarded bam
# sort and then move/rename BAM
samtools sort input_files/discarded.bam  sorted.discarded 
mv sorted.discarded.bam out/discarded_BAM/output/$samplename.refined.primerclippeddiscarded.bam 
# index sorted bam and move index into own output folder
samtools index out/discarded_BAM/output/$samplename.refined.primerclippeddiscarded.bam 
mv out/discarded_BAM/output/$samplename.refined.primerclippeddiscarded.bam.bai out/discarded_BAM_BAI/output/$samplename.refined.primerclippeddiscarded.bam.bai

# move the metrics file into the output folder
mv input_files/$samplename.refined.primerclipped.metrics out/metrics/QC/

# upload all outputs
dx-upload-all-outputs --parallel
