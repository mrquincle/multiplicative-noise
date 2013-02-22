<!-- Uses markdown syntax for neat display at github -->

# Multiplicative noise

Self-organized criticality (SOC) is an interesting phenomenon that occurs in sandpiles. It may occur on other systems too such as earthquake aftershock distributions. The hallmark of SOC is avalanches, which however can also be caused by a wide variety of other phenomena that have nothing to do with SOC. Ivan et al. describe Langevin equations for sandpiles in the case of grains dissipating (getting lost) in the bulk versus only leaving the system at the sides. In the latter case there is a nice power-law distribution, in the former case, there is no such thing. The Langevin equations do indeed demonstrate this. To solve those Langevin equations with multiplicative noise, an advanced splitting scheme is used by the authors. I would like to thank Ivan for providing me the Fortran code for that, it made writing it in C++ much easier.

This code did run with a former Boost version flawlessly in the multi-thread setting. However, the current Boost version seems to choke on it.

## Copyrights
The copyrights (2012) belong to:

- Author: Anne van Rossum
- Author: Ivan Dornic
- Almende B.V., http://www.almende.com and DO bots B.V., http://www.dobots.nl
- Rotterdam, The Netherlands


