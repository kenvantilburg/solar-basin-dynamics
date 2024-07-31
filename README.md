# solar-basin-dynamics

Code repository for the paper **Orbital Dynamics of Solar Basin Particles** by Cara Giovanetti, Robert Lasenby, and Ken Van Tilburg.

## Abstract

We study the dynamics of the solar basin---the accumulated population of weakly-interacting particles on bound orbits in the Solar System.
We focus on particles starting off on Sun-crossing orbits, corresponding to initial conditions of production inside the Sun, and investigate their evolution over the age of the Solar System.
A combination of analytic methods, secular perturbation theory, and direct numerical integration of orbits sheds light on the long- and short-term evolution of a population of test particles orbiting the Sun and perturbed by the planets.
Our main results are that the effective lifetime of a solar basin at Earth's location is $\tau_{\rm eff} = 1.20\pm 0.09 \,\mathrm{Gyr}$, and that there is annual (semi-annual) modulation of the basin density with known phase and amplitude at the fractional level of 6.5\% (2.2\%). 
These results have important implications for direct detection searches of solar basin particles, and the strong temporal modulation signature yields a robust discovery channel.
Our simulations can also be interpreted in the context of gravitational capture of dark matter in the Solar System, with consequences for any dark-matter phenomenon that may occur below the local escape velocity.

## Code

The dependencies of the code are listed in [environments.yml](environment.yml).

The [code](code/) folder contains various Jupyter and Mathematica notebooks that reproduce the plots in the paper. These are linked to from the paper.

## Authors

-  Cara Giovanetti (cg3566@nyu.edu)
-  Robert Lasenby (rlasenby@stanford.edu)
-  Ken Van Tilburg (kenvt@nyu.edu)

## Citation

If you use this code, please cite our paper:
```
[put Bibtex here]
```
and you might want to cite the original "stellar basins" papers for axions and dark photons, respectively:
```
@article{VanTilburg:2020jvl,
    author = "Van Tilburg, Ken",
    title = "{Stellar basins of gravitationally bound particles}",
    eprint = "2006.12431",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1103/PhysRevD.104.023019",
    journal = "Phys. Rev. D",
    volume = "104",
    number = "2",
    pages = "023019",
    year = "2021"
}
```
```
@article{Lasenby:2020goo,
    author = "Lasenby, Robert and Van Tilburg, Ken",
    title = "{Dark photons in the solar basin}",
    eprint = "2008.08594",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    doi = "10.1103/PhysRevD.104.023020",
    journal = "Phys. Rev. D",
    volume = "104",
    number = "2",
    pages = "023020",
    year = "2021"
}
```
