# zk-Stark

This repository is created to learn implementations of zero-knowledge
STARK protocol.

## Credits

Foundation of this code relies on [great tutorial by aszepieniec](https://aszepieniec.github.io/stark-anatomy/).

## TODO

- [ ] There is no compression of transition constraints into one master constraint
- [ ] In production systems the length of the codeword is often reduced not by a factor 2 but a small power of 2. This optimization reduces the proof size and might even generate running time improvements
- [ ] Arithmetization step, which transforms the initial computation into an arithmetic constraint system
- [ ] Extension fields. [EthSTARK](https://github.com/starkware-libs/ethSTARK) operates over the finite field defined by a 62-bit prime, but the FRI step operates over a quadratic extension field thereof in order to target a higher security level.
- [ ] Advanced STARKs may define more constraint types (referring to AIR) in order to deal with memory or with consistency of registers within one cycle
