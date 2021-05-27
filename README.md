# Artificial PDB coordinates
Let's make non-existing PDB coordinates which describe any geometric shape and for any taste.

## PDB file format

<pre class="bold fixed">column            1         2         3         4         5         6         7         8
number   12345678901234567890123456789012345678901234567890123456789012345678901234567890

         <span class="f1">ATOM  </span><span class="f2">    1</span> <span class="f3"> N  </span><span class="f4"> </span><span class="f5">GLY</span> <span class="f6">A</span><span class="f7">   3</span><span class="f8"> </span>   <span class="f9">  17.119</span><span class="f10">   0.186</span><span class="f11">  36.320</span><span class="f12">  1.00</span><span class="f13"> 64.10</span>          <span class="f14"> N</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    2</span> <span class="f3"> CA </span><span class="f4"> </span><span class="f5">GLY</span> <span class="f6">A</span><span class="f7">   3</span><span class="f8"> </span>   <span class="f9">  16.944</span><span class="f10">  -0.800</span><span class="f11">  35.208</span><span class="f12">  1.00</span><span class="f13"> 63.46</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    3</span> <span class="f3"> C  </span><span class="f4"> </span><span class="f5">GLY</span> <span class="f6">A</span><span class="f7">   3</span><span class="f8"> </span>   <span class="f9">  16.818</span><span class="f10">  -0.087</span><span class="f11">  33.851</span><span class="f12">  1.00</span><span class="f13"> 61.22</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    4</span> <span class="f3"> O  </span><span class="f4"> </span><span class="f5">GLY</span> <span class="f6">A</span><span class="f7">   3</span><span class="f8"> </span>   <span class="f9">  15.721</span><span class="f10">   0.337</span><span class="f11">  33.463</span><span class="f12">  1.00</span><span class="f13"> 62.81</span>          <span class="f14"> O</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    5</span> <span class="f3"> N  </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  17.944</span><span class="f10">   0.077</span><span class="f11">  33.129</span><span class="f12">  1.00</span><span class="f13"> 57.39</span>          <span class="f14"> N</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    6</span> <span class="f3"> CA </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  17.950</span><span class="f10">   0.742</span><span class="f11">  31.815</span><span class="f12">  1.00</span><span class="f13"> 53.27</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    7</span> <span class="f3"> C  </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  18.005</span><span class="f10">  -0.247</span><span class="f11">  30.629</span><span class="f12">  1.00</span><span class="f13"> 49.78</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    8</span> <span class="f3"> O  </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  19.086</span><span class="f10">  -0.678</span><span class="f11">  30.218</span><span class="f12">  1.00</span><span class="f13"> 48.17</span>          <span class="f14"> O</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">    9</span> <span class="f3"> CB </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  19.191</span><span class="f10">   1.613</span><span class="f11">  31.898</span><span class="f12">  1.00</span><span class="f13"> 54.33</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">   10</span> <span class="f3"> CG </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  20.161</span><span class="f10">   0.686</span><span class="f11">  32.625</span><span class="f12">  1.00</span><span class="f13"> 55.45</span>          <span class="f14"> C</span><span class="f15">  </span>
         <span class="f1">ATOM  </span><span class="f2">   11</span> <span class="f3"> CD </span><span class="f4"> </span><span class="f5">PRO</span> <span class="f6">A</span><span class="f7">   4</span><span class="f8"> </span>   <span class="f9">  19.305</span><span class="f10">   0.019</span><span class="f11">  33.701</span><span class="f12">  1.00</span><span class="f13"> 55.83</span>          <span class="f14"> C</span><span class="f15">  </span>

field id <span class="f1">  1   </span><span class="f2">  2  </span> <span class="f3">  3 </span><span class="f4">4</span><span class="f5"> 5 </span> <span class="f6">6</span><span class="f7">  7 </span><span class="f8">8</span>   <span class="f9">    9   </span><span class="f10">   10   </span><span class="f11">   11   </span><span class="f12">  12  </span><span class="f13">  13  </span>          <span class="f14">14</span><span class="f15">15</span>
</pre>
<table class="centered">
<tbody><tr>
    <td class="header">field id</td>
    <td class="header left">definition</td>
    <td class="header">length</td>
    <td class="header">format</td>
    <td class="header">range</td>
    <td class="header left">string slicing (Python)</td>
</tr>

<tr>
    <td class="f1 fixed bold">1</td>
    <td class="left"><tt>"ATOM  "</tt> or <tt>"HETATM"</tt></td>
    <td>6</td>
    <td class="fixed">{:6s}</td>
    <td>01-06</td>
    <td class="left">[0:6]</td>
</tr>
<tr>
    <td class="f2 fixed bold">2</td>
    <td class="left">atom serial number</td>
    <td>5</td>
    <td class="fixed">{:5d}</td>
    <td>07-11</td>
    <td class="left">[6:11]</td>
</tr>
<tr>
    <td class="f3 fixed bold">3</td>
    <td class="left">atom name</td>
    <td>4</td>
    <td class="fixed">{:^4s}</td>
    <td>13-16</td>
    <td class="left">[12:16]</td>
</tr>
<tr>
    <td class="f4 fixed bold">4</td>
    <td class="left">alternate location indicator</td>
    <td>1</td>
    <td class="fixed">{:1s}</td>
    <td>17</td>
    <td class="left">[16:17]</td>
</tr>
<tr>
    <td class="f5 fixed bold">5</td>
    <td class="left">residue name</td>
    <td>3</td>
    <td class="fixed">{:3s}</td>
    <td>18-20</td>
    <td class="left">[17:20]</td>
</tr>
<tr>
    <td class="f6 fixed bold">6</td>
    <td class="left">chain identifier</td>
    <td>1</td>
    <td class="fixed">{:1s}</td>
    <td>22</td>
    <td class="left">[21:22]</td>
</tr>
<tr>
    <td class="f7 fixed bold">7</td>
    <td class="left">residue sequence number</td>
    <td>4</td>
    <td class="fixed">{:4d}</td>
    <td>23-26</td>
    <td class="left">[22:26]</td>
</tr>
<tr>
    <td class="f8 fixed bold">8</td>
    <td class="left">code for insertion of residues</td>
    <td>1</td>
    <td class="fixed">{:1s}</td>
    <td>27</td>
    <td class="left">[26:27]</td>
</tr>
<tr>
    <td class="f9 fixed bold">9</td>
    <td class="left">orthogonal coordinates for X (in Angstroms)</td>
    <td>8</td>
    <td class="fixed">{:8.3f}</td>
    <td>31-38</td>
    <td class="left">[30:38]</td>
</tr>
<tr>
    <td class="f10 fixed bold">10</td>
    <td class="left">orthogonal coordinates for Y (in Angstroms)</td>
    <td>8</td>
    <td class="fixed">{:8.3f}</td>
    <td>39-46</td>
    <td class="left">[38:46]</td>
</tr>
<tr>
    <td class="f11 fixed bold">11</td>
    <td class="left">orthogonal coordinates for Z (in Angstroms)</td>
    <td>8</td>
    <td class="fixed">{:8.3f}</td>
    <td>47-54</td>
    <td class="left">[46:54]</td>
</tr>
<tr>
    <td class="f12 fixed bold">12</td>
    <td class="left">occupancy</td>
    <td>6</td>
    <td class="fixed">{:6.2f}</td>
    <td>55-60</td>
    <td class="left">[54:60]</td>
</tr>
<tr>
    <td class="f13 fixed bold">13</td>
    <td class="left">temperature factor</td>
    <td>6</td>
    <td class="fixed">{:6.2f}</td>
    <td>61-66</td>
    <td class="left">[60:66]</td>
</tr>
<tr>
    <td class="f14 fixed bold">14</td>
    <td class="left">element symbol</td>
    <td>2</td>
    <td class="fixed">{:&gt;2s}</td>
    <td>77-78</td>
    <td class="left">[76:78]</td>
</tr>
<tr>
    <td class="f15 fixed bold">15</td>
    <td class="left">charge on the atom</td>
    <td>2</td>
    <td class="fixed">{:2s}</td>
    <td>79-80</td>
    <td class="left">[78:80]</td>
</tr>
</tbody></table>

`ATOM` line in PDB file must look like:
```python
"{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)
```

- For all protein models use `ATOM`, but `HETATM` for ligands.
- Atom name = `CA`
- Residue name for `ATOM` = `ALA`, but for `HETATM` can vary, e.g. `LIG`, `LG1`, `LG2`.
- Occupancy = 1.00
- Temperature factor = 0.00
- Element symbol = `C`
- Charge on the atom = absent
## References

- [PDB file format](https://cupnet.net/pdb-format/)
- [PDB to FASTA in Bash](https://cupnet.net/pdb2fasta/)

Project Organization
------------

    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── config.yaml        <- A configuration file for the scripts.
    |
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    |
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    |
    ├── ccenv.yml          <- A Conda environmet file.
    |
    ├── misc
    |
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    |
    ├── output
    |
    ├── pdb
    │   ├── processed
    │   └── raw
    |
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    |
    ├── src
    │   ├── data           <- Scripts to download or generate data
    │   │   ├── make_dataset.py
    │   │   └── make_pdb.py
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    |
    ├── tmp                <- for temporary files
    |
    ├── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io
    |
    └── .env               <- Stores environment variables for Dotenv module. Do not track with version control.

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
