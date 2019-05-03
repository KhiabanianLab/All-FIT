
# All-FIT

All-FIT - *Allele-Frequency-based Imputation of Tumor Purity* infers specimen purity from tumor-only samples sequenced with deep sequencing. It is developed in [Khiabanian Lab](http://khiabanian-lab.org/) by Jui Wan Loh and Hossein Khiabanian.<br/>
All-FIT is implemented in Python 3.

## How to run

### Dependencies

The following Python package is required:

- [numpy](http://www.numpy.org)
- [matplotlib](https://matplotlib.org/)
- [scipy](https://www.scipy.org/install.html)
- [seaborn](https://pypi.org/project/seaborn/)

### Installation

Just clone this repository and run All-FIT as described below.

### Usage example

Running All-FIT on test/input/sampleFile1.txt with columns containing uniqueID, variant allele frequencies, depth and ploidy information:

```
python All-FIT.py -i All-FIT/test/input/sampleFile1.txt -d All-FIT/test/output -o sampleFile1
sampleFile1     0.78    0.77,0.78,0.79
```

Inferred purity is printed on stardard output as shown above (estimated purity = 0.78 with two standard deviation ranges from 0.77 to 0.79); one text file and two png files are generated in All-FIT/test/output.<br/>

Try:

```
python All-FIT.py  --help
```

for detailed description of input parameters


