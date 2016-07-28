/home/tbecker/software/varsim/varsim.py \
--vc_in_vcf /data/tbecker/varsim_data/All.vcf.gz \
--sv_insert_seq /data/tbecker/varsim_data/insert_seq.txt \
--sv_dgv /data/tbecker/varsim_data/GRCh37_hg19_supportingvariants_2013-07-23.txt \
--reference /data/tbecker/varsim_data/hs37d5.fa --id simu \
--read_length 100 --mean_fragment_size 350 --sd_fragment_size 50 --nlanes 5 --total_coverage 1 \
--vc_num_snp 3000000 --vc_num_ins 100000 --vc_num_del 100000 --vc_num_mnp 50000 --vc_num_complex 50000 \
--vc_percent_novel 0.01 --vc_min_length_lim 0 --vc_max_length_lim 49 \
--sv_num_ins 2000 --sv_num_del 2000 --sv_num_dup 200 --sv_num_inv 1000 --sv_percent_novel 0.01 \
--sv_min_length_lim 50 --sv_max_length_lim 1000000 \
--simulator_executable /home/tbecker/software/art_bin_VanillaIceCream/art_illumina --simulator art \
--out_dir /data/tbecker/varsim_test --log_dir /data/tbecker/varsim_test --work_dir /data/tbecker/varsim_test
#--vcfs [Optional VCF file to include, variants.vcf]