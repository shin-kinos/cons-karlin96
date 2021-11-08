# cons-karlin96 

A Rust program that scores residue conservation a site in a Multiple Sequence Alignment (MSA) by using Sum-of-pairs measure (SP measure). 

# Description 

* The program scores residue conservation in a MSA. 
* The scoring method is based on SP measure. 
* It takes account of normalized substitution matrices and gap penalty. 

## Dependencies 

* `colored` ( https://github.com/mackwic/colored ) 

``` 
[dependencies]
colored = "2.0"
``` 

## Installation 

You can compile this program by using `Cargo`. 📦🦀

[e. g.] 

``` 
% cd cons-karlin96-main
% cargo build --release
``` 
Then the object file is generated in `./target/release` directory. 

## Scoring method 

### Conservation score 
The conservation score is calculated based on SP measure as follows [1]: 

<img width="1440" alt="cons-karlin96_equation01" src="https://user-images.githubusercontent.com/83740080/140666320-7d079490-5ad3-4d20-af18-3de74bb3b1dd.png"> 

where ***N*** is the length of a site (or number of the sequences in the MSA ), ***Sx( i )*** is an amino acid at site ***i*** in sequence ***x*** and ***M*** is a normalized substitution matrix. 

The normalized substitution ***M*** is given based on Karlin-like method as follows: 

<img width="1440" alt="cons-karlin96_equation02" src="https://user-images.githubusercontent.com/83740080/140666340-8e856753-9ceb-4368-95e8-d3580c86aa6c.png"> 

where ***m*** is a typical substitution scoring matrix such as BLOSUM series. This normalization method modify a matrix so that all the diagonal elements are the maxima at 1. If ***a*** or ***b*** is gap, ***M(a, b)*** gets the minimum element of the normalized matrix as gap penalty. 

### Substitution scoring matrices 

This program supports 11 substitution scoring matrices as follows [1]:

* BLOSUM45 [2]
* BLOSUM50 
* BLOSUM62 
* BLOSUM70 
* BLOSUM80 
* BLOSUM90 
* PAM30 
* PAM70 
* PAM250 [3] 
* Modified version of PET91* [4]
* Modified version of BLOSUM62*

*These matrices are modified so that all of their diagonal elements are constant at the maximum (`W vs W = 11` in BLOSUM62 and `W vs W = 15` in PET91). 

### ⚠️ Sequence weighting 

This program does **NOT** support sequence weighting ! 

## Input file format 

Aligned Multi-FASTA format. NOTE that nucleotide sequences are not supported. 

See some example input files in `demo` directory. 

## Usage 

Major arguments:

`-i` : Input filename in aligned Multi-FASTA format, REQUIRED.

`-o` : Output filename, REQUIRED.

`-m` : Substitution scoring matrices (default "blosum62"). 

[e. g.]

``` 
% ./cons-karlin96 -i input.fasta -o output.txt -m pet91mod -c yes -t no
``` 

Type `-h` to see other available options. 

## Output 

Site number `\t` Conservation score `\t` Composition of the site

[e.g.] 

<img width="946" alt="cons-karlin96_output" src="https://user-images.githubusercontent.com/83740080/140674282-0cd5f8ee-ba9f-4785-8c3d-1ea7b5b49e44.png"> 

## References 

1. Karlin, Samuel, and Luciano Brocchieri. "Evolutionary conservation of RecA genes in relation to protein structure and function." Journal of bacteriology 178.7 (1996): 1881-1894. 
2. Henikoff, Steven, and Jorja G. Henikoff. "Amino acid substitution matrices from protein blocks." Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919. 
3. Dayhoff, M., R. Schwartz, and B. Orcutt. "A model of evolutionary change in proteins." Atlas of protein sequence and structure 5 (1978): 345-352. 
4. Jones, David T., William R. Taylor, and Janet M. Thornton. "The rapid generation of mutation data matrices from protein sequences." Bioinformatics 8.3 (1992): 275-282. 

