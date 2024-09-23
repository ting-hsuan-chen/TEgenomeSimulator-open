A summary of the genome simulated by TEgenomeSimulator
================
Ting-Hsuan Chen
05 March, 2024

# Simulation information

**Simulated genome: donghong**
**Range of simulated TE loci per TE family: 5-10**

|  CPU| RAM       | Time     | Timestamp\_of\_completion |
|----:|:----------|:---------|:--------------------------|
|    1| 842.30 MB | 00:00:51 | 2024-02-21 14:14:22       |

In this simulation, TEgenomeSimulator did the followings:

1.  It took the TE library file, **combined\_curated\_TE\_lib\_ATOSZM\_selected.fasta**, which comprises curated TE family sequences from ***A. thalian***, ***Z. maize*** and ***O. sativa***, and simulates multiple TE copies with sequence variations depending on the parameters specified in the table: **TElib\_sim\_list\_5\_10.table**. In addition to nucleotide substitution and INDEL, it also simulated fragmentation (as a proportion of TE truncated from 5' end), nested insertion (only for Copia and Gypsy), as well as target site duplication.

2.  The simulated TE copies were then randomly inserted into the user-provided TE-depleted genome, **Donghong.chromosomes.only.fa.nonTE**, where TE sequences had been exhaustively detected by multiple TE annotators (e.g. EDTA, RepeatModeler, and EarlGrey) and removed. The final simulated genome can be utilised for benchmarking TE annotators.

# Loaded files for creating this report

-   Genome fasta index file: donghong\_5\_10\_genome\_sequence\_out\_nest.fasta.fai
-   All TE fasta index file: donghong\_5\_10\_repeat\_sequence\_out\_nest.fasta.fai
-   All TE gff file: donghong\_5\_10\_repeat\_annotation\_out\_nest.gff

# The proportion of the TE/nonTE sequence

-   Total simulated genome size: 306,230,938 bp
-   Total simulated TE bases: 22,638,471 bp

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig1-1.png" width="50%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 001\_sim\_genome\_size\_donghong\_5\_10.png

# Breaking down the simulated TEs by superfamily

There are total 16 TE superfamilies included in the curated TE library. The following table shows the number of loci, bases and family of each superfamily, as well as the percentage of loci, bases and family.

<table>
<colgroup>
<col width="18%" />
<col width="7%" />
<col width="10%" />
<col width="14%" />
<col width="17%" />
<col width="15%" />
<col width="19%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">TE_superfamily</th>
<th align="right">loci</th>
<th align="right">bp</th>
<th align="right">family_count</th>
<th align="right">loci_percentage</th>
<th align="right">bp_percentage</th>
<th align="right">family_percentage</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">01_LTR/Copia</td>
<td align="right">920</td>
<td align="right">2907611</td>
<td align="right">139</td>
<td align="right">7.51</td>
<td align="right">12.84</td>
<td align="right">8.44</td>
</tr>
<tr class="even">
<td align="left">02_LTR/Gypsy</td>
<td align="right">753</td>
<td align="right">2476491</td>
<td align="right">109</td>
<td align="right">6.15</td>
<td align="right">10.94</td>
<td align="right">6.62</td>
</tr>
<tr class="odd">
<td align="left">03_LTR/Solo</td>
<td align="right">21</td>
<td align="right">28080</td>
<td align="right">3</td>
<td align="right">0.17</td>
<td align="right">0.12</td>
<td align="right">0.18</td>
</tr>
<tr class="even">
<td align="left">04_LTR/unknown</td>
<td align="right">84</td>
<td align="right">110541</td>
<td align="right">11</td>
<td align="right">0.69</td>
<td align="right">0.49</td>
<td align="right">0.67</td>
</tr>
<tr class="odd">
<td align="left">05_LINE/L1</td>
<td align="right">89</td>
<td align="right">312729</td>
<td align="right">12</td>
<td align="right">0.73</td>
<td align="right">1.38</td>
<td align="right">0.73</td>
</tr>
<tr class="even">
<td align="left">06_LINE/unknown</td>
<td align="right">710</td>
<td align="right">1470449</td>
<td align="right">94</td>
<td align="right">5.80</td>
<td align="right">6.50</td>
<td align="right">5.71</td>
</tr>
<tr class="odd">
<td align="left">07_SINE/tRNA</td>
<td align="right">38</td>
<td align="right">8237</td>
<td align="right">5</td>
<td align="right">0.31</td>
<td align="right">0.04</td>
<td align="right">0.30</td>
</tr>
<tr class="even">
<td align="left">08_SINE/unknown</td>
<td align="right">307</td>
<td align="right">153493</td>
<td align="right">43</td>
<td align="right">2.51</td>
<td align="right">0.68</td>
<td align="right">2.61</td>
</tr>
<tr class="odd">
<td align="left">09_DNA/CACTA</td>
<td align="right">789</td>
<td align="right">2109253</td>
<td align="right">110</td>
<td align="right">6.44</td>
<td align="right">9.32</td>
<td align="right">6.68</td>
</tr>
<tr class="even">
<td align="left">10_DNA/hAT</td>
<td align="right">1891</td>
<td align="right">1627563</td>
<td align="right">249</td>
<td align="right">15.44</td>
<td align="right">7.19</td>
<td align="right">15.12</td>
</tr>
<tr class="odd">
<td align="left">11_DNA/MuDR</td>
<td align="right">2993</td>
<td align="right">6228941</td>
<td align="right">393</td>
<td align="right">24.43</td>
<td align="right">27.51</td>
<td align="right">23.86</td>
</tr>
<tr class="even">
<td align="left">12_DNA/Harbinger</td>
<td align="right">23</td>
<td align="right">68688</td>
<td align="right">3</td>
<td align="right">0.19</td>
<td align="right">0.30</td>
<td align="right">0.18</td>
</tr>
<tr class="odd">
<td align="left">13_DNA/Mariner</td>
<td align="right">58</td>
<td align="right">32961</td>
<td align="right">7</td>
<td align="right">0.47</td>
<td align="right">0.15</td>
<td align="right">0.43</td>
</tr>
<tr class="even">
<td align="left">14_RC/Helitron</td>
<td align="right">2424</td>
<td align="right">4835648</td>
<td align="right">314</td>
<td align="right">19.79</td>
<td align="right">21.36</td>
<td align="right">19.06</td>
</tr>
<tr class="odd">
<td align="left">15_MITE/Stow</td>
<td align="right">286</td>
<td align="right">60556</td>
<td align="right">40</td>
<td align="right">2.33</td>
<td align="right">0.27</td>
<td align="right">2.43</td>
</tr>
<tr class="even">
<td align="left">16_MITE/Tourist</td>
<td align="right">865</td>
<td align="right">207940</td>
<td align="right">115</td>
<td align="right">7.06</td>
<td align="right">0.92</td>
<td align="right">6.98</td>
</tr>
<tr class="odd">
<td align="left">Total</td>
<td align="right">12251</td>
<td align="right">22639181</td>
<td align="right">1647</td>
<td align="right">100.00</td>
<td align="right">100.00</td>
<td align="right">100.00</td>
</tr>
</tbody>
</table>

Table saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: sim\_genome\_summary\_TEsuperfamily\_donghong\_5\_10.csv

## Total simulated TE loci categorised by TE superfamily

-   Total simulated TE loci: 12,251 loci

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig2-1.png" width="60%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 002\_sim\_te\_loci\_donghong\_5\_10.png

## Total simulated TE bases categorised by TE superfamily

-   Total simulated TE bases: 22,639,181 bp

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig3-1.png" width="60%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 003\_sim\_te\_bp\_donghong\_5\_10.png

## Simulated TE family per superfamily

This part only depends on the curated TE library. It shouldn't make any difference between simulations.

-   Total simulated TE families: 1647 families

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig4-1.png" width="60%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 004\_sim\_te\_family\_donghong\_5\_10.png

# Extracting full information from TE's gff file

Have a look at the full info extracted from TE's gff file (row 100 to 109):

Full table saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: sim\_genome\_TE\_insertion\_info\_donghong\_5\_10.csv

# Simulation of nested TE insertions

## Nested and non-nested Copia

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig5-1.png" width="50%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 005\_sim\_te\_nested\_copia\_donghong\_5\_10.png

## TE loci cut by nested Copia

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig6-1.png" width="80%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 006\_sim\_te\_cutby\_nested\_copia\_donghong\_5\_10.png

## Nested and non-nested Gypsy

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig7-1.png" width="50%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 007\_sim\_te\_nested\_gypsy\_donghong\_5\_10.png

## TE loci cut by nested Gypsy

<img src="09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig8-1.png" width="80%" />

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 008\_sim\_te\_cutby\_nested\_gypsy\_donghong\_5\_10.png

# The distribution of TE loci identity

How did TEgenomeSimulator simulate sequence identity in this simulation?

1.  The distribution of sequence identity of each TE family was defined by the **idn** and **sd** values (user provided) stored in **TElib\_sim\_list\_5\_10.table**.

2.  TEgenomeSimulator took the **idn** as mean identity and **sd** as standard deviation to create a distrubution, from which a value was sampled as the **simulated identity** of a TE member. Therefore, **simulated divergence = 1 - simulated identity**

3.  Inherited from it's predecessor, [denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval), TEgenomSimulator broke down the **simulated divergence** into **substitution** and **INDELs** (i.e. **divergence = substitution % + INDEL %**).

4.  The **INDEL %** was defined by the **indels** value from **TElib\_sim\_list\_5\_10.table**.

![](09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig9-1.png)

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 009\_sim\_te\_loci\_identity\_donghong\_5\_10.png

# The distribution of TE loci integrity

How did TEgenomeSimulator simulate sequence integrity in this simulation?

1.  TEgenomeSimulator considers sequence integrity as **integrity = (1 - (TE locus length / full length))**

2.  The length of a TE locus is decided by INDELs and fragmentation.

3.  In the fragmentation step, if a TE length is shorter than 500 bp, TEgenomeSimulator would randomly select a value between 70 and 90 as the fraction of the sequence to be removed from 5' end; otherwise the simulator randomly chooses a value between 40 and 90. (This is the same as in [denovoTE-eval](https://github.com/IOB-Muenster/denovoTE-eval))

4.  The number of fragmented loci was defined by the value of **frag** in **TElib\_sim\_list\_5\_10.table**. This value was taken as a proportion of total TE loci to undergo fragmentation step.

![](09_summarise_sim_genome_DH_5_10_files/figure-markdown_github/fig10-1.png)

High quality image saved in
/workspace/cflthc/scratch/2022\_Actinidia\_TE/10.01\_TEgenomeSimulator\_output/sim\_donghong\_5\_10\_genome/report/

File name: 010\_sim\_te\_loci\_integrity\_donghong\_5\_10.png
