from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Gide-data/main/"

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"))
    resources:
        mem_mb=3000,
        disk_mb=3000
    shell:
        """
        Rscript -e \
        '
        load(paste0("{prefix}", "annotation/Gencode.v40.annotation.RData"))
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "Gide", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        """

rule format_data:
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv")
    input:
        S3.remote(prefix + "download/mel_gide19_exp_data.csv"),
        S3.remote(prefix + "download/mel_gide19_survival_data.csv"),
        S3.remote(prefix + "download/mel_gide19_clin_data.csv"),
        S3.remote(prefix + "download/CLIN_GIDE.txt")
    shell:
        """
        Rscript scripts/Format_Data.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_data:
    output:
        S3.remote(prefix + "download/mel_gide19_exp_data.csv"),
        S3.remote(prefix + "download/mel_gide19_survival_data.csv"),
        S3.remote(prefix + "download/mel_gide19_clin_data.csv"),
        S3.remote(prefix + "download/CLIN_GIDE.txt")
        
    shell:
        """
        wget {data_source}mel_gide19_exp_data.csv -O {prefix}download/mel_gide19_exp_data.csv
        wget -O {prefix}download/mel_gide19_survival_data.csv {data_source}mel_gide19_survival_data.csv
        wget -O {prefix}download/mel_gide19_clin_data.csv {data_source}mel_gide19_clin_data.csv
        wget -O {prefix}download/CLIN_GIDE.txt {data_source}CLIN_GIDE.txt
        """ 