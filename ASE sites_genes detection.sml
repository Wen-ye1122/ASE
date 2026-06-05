sample="" # sample id
# 1. ASE sites detection using phaser (python2 is needed)
rule ASE detection:
  input:
      phASER_dir="", # Path to STAR,
      vcf="{sample}"
  output:
        ""
    threads: 5
    shell:
      ```
      if [ pair end ]  # I don't know how to set or define ???
      ```
