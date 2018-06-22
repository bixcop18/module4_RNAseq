#!/bin/sh
#$ -S /bin/bash
#$ -v PATH=/home/oasis/data/webcomp/RAMMCAP-ann/bin:/usr/local/bin:/usr/bin:/bin
#$ -v BLASTMAT=/home/oasis/data/webcomp/RAMMCAP-ann/blast/bin/data
#$ -v LD_LIBRARY_PATH=/home/oasis/data/webcomp/RAMMCAP-ann/gnuplot-install/lib
#$ -v PERL5LIB=/home/hying/programs/Perl_Lib
#$ -q all.q
#$ -pe orte 2


#$ -e /home/oasis/data/webcomp/web-session/1529553255/1529553255.err
#$ -o /home/oasis/data/webcomp/web-session/1529553255/1529553255.out
cd /home/oasis/data/webcomp/web-session/1529553255
faa_stat.pl 1529553255.fas.0

/home/oasis/data/NGS-ann-project/apps/cd-hit/cd-hit-est -i 1529553255.fas.0 -d 0 -o 1529553255.fas.1 -c 0.98 -n 10 -l 11  -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000
faa_stat.pl 1529553255.fas.1
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1529553255.fas.1.clstr > 1529553255.fas.1.clstr.sorted
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list.pl 1529553255.fas.1.clstr 1529553255.clstr.dump
gnuplot1.pl < 1529553255.fas.1.clstr > 1529553255.fas.1.clstr.1; gnuplot2.pl 1529553255.fas.1.clstr.1 1529553255.fas.1.clstr.1.png
/home/oasis/data/NGS-ann-project/apps/cd-hit/cd-hit-est -i 1529553255.fas.1 -d 0 -o 1529553255.fas.2 -c 0.95 -n 10 -l 11  -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000
faa_stat.pl 1529553255.fas.2
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1529553255.fas.2.clstr > 1529553255.fas.2.clstr.sorted
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list.pl 1529553255.fas.2.clstr 1529553255.clstr.dump
gnuplot1.pl < 1529553255.fas.2.clstr > 1529553255.fas.2.clstr.1; gnuplot2.pl 1529553255.fas.2.clstr.1 1529553255.fas.2.clstr.1.png
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_rev.pl 1529553255.fas.1.clstr 1529553255.fas.2.clstr > 1529553255.fas.2-0.clstr
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1529553255.fas.2-0.clstr > 1529553255.fas.2-0.clstr.sorted
/home/oasis/data/NGS-ann-project/apps/cd-hit/cd-hit-est -i 1529553255.fas.2 -d 0 -o 1529553255.fas.3 -c 0.9 -n 8  -r 1 -G 1 -g 1 -b 20 -s 0.0 -aL 0.0 -aS 0.0 -T 4 -M 32000
faa_stat.pl 1529553255.fas.3
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1529553255.fas.3.clstr > 1529553255.fas.3.clstr.sorted
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list.pl 1529553255.fas.3.clstr 1529553255.clstr.dump
gnuplot1.pl < 1529553255.fas.3.clstr > 1529553255.fas.3.clstr.1; gnuplot2.pl 1529553255.fas.3.clstr.1 1529553255.fas.3.clstr.1.png
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_rev.pl 1529553255.fas.2-0.clstr 1529553255.fas.3.clstr > 1529553255.fas.3-0.clstr
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_sort_by.pl no < 1529553255.fas.3-0.clstr > 1529553255.fas.3-0.clstr.sorted
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1529553255.clstr.dump 1529553255.clstr_no.dump
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1529553255.clstr.dump 1529553255.clstr_len.dump len
/home/oasis/data/NGS-ann-project/apps/cd-hit/clstr_list_sort.pl 1529553255.clstr.dump 1529553255.clstr_des.dump des
tar -zcf 1529553255.result.tar.gz * --exclude=*.dump --exclude=*.env
echo hello > 1529553255.ok
