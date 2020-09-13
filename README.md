# Mock Fungible genetic simulation code

The easy way to fuse two attributes.

# API format

```
GET https://cryptopepes.io/api/getPepe/1234
```

(Replace 1234 with MFRG chip ID)

Retrieved:

```
{
  "gen": 0,
  "genotype": "8d58f1cd1858b473bbe38697eef629cb292e592bd2ca9801748b76e3cdc88707a9f2af5e826bd695c3611d3ca08fb5c3a6a13c67009412fc36b39865cba76567",
  "mother": "0",
  "father": "0",
  "master": "0x831491b5cf7ba0d626672abf71a1420f7d3844ac",
  ... *unimportant data*
  "pepeId": "1"
}
```

## Genotype format

```
2 strings, each string is 256 bits
Each string is 2 times 128 bits (a chromosome)
Each chromosome consists of genes.
Each gene has a length, and a starting position (read left to right)

|                           string 0                           ||                           string 1                           |
|         chromosome 0 (A)     ||        chromosome 1  (B)     ||         chromosome 0   (C)   ||           chromosome 1  (D)  |
8d58f1cd1858b473bbe38697eef629cb292e592bd2ca9801748b76e3cdc88707a9f2af5e826bd695c3611d3ca08fb5c3a6a13c67009412fc36b39865cba76567

Now split, and re-order it to:

string P: BA
string Q: DC

For Web3 users:
If you're using the output of the getMFRG() smartcontract call, then dna[0] == P, and dna[1] == Q
```


`P` and `Q` have the same interpretation mechanic, but there's one difference: `Q` is CAN be DOMINANT over `P`.
This means that if an allel on `Q` is dominant, it will be used for the expression. If it's not, then `P` will be the expressed DNA.


## Expression estimation

`output.json` defines all the *counts* (not direct chances) of each property showing up for a certain bit being `0` (`c0`) or `1` (`c1`) on `P`.
It also does the same for the other side `Q`, but instead it's `c0dom` and `c1dom`. There's a special property `X` that sums all the times that a `1`/`0` on `Q` was not dominant.

### Math

#### Dependency check

The chance of a certain property `T` showing up, when we know it's already there. This can be used to check how independent the bits are.

In pseudo code:
```
spec = outputJSON array

totalChance = 1

for each DNA bit "i" [range [0...255]] (
    var chanceP;
    if currentDNA_P[i] == 1 (
        chanceP = spec[i]["c1"][my_prop] || 0     <--- 0 if my_prop key is absent
    ) else (
        chanceP = spec[i]["c0"][my_prop] || 0
    )
    var chanceQ;
    if currentDNA_Q[i] == 1 (
        chanceQ = spec[i]["c1dom"][my_prop] || 0
    ) else (
        chanceQ = spec[i]["c0dom"][my_prop] || 0
    )
    var maxCount = spec[i]["max"]
    # Normalize the counts, now they're chances, range [0...1]
    chanceP = chanceP / maxCount
    chanceQ = chanceQ / maxCount
    
    # This is just the chance for this single bit contributing to the property being expressed
    var chance = chanceQ + chanceP
    
    // If there's no chance for either 0 or 1, recessive (string P),
    //  then this bit does not affect the expression of the property.
    if spec[i]["c0"].contains(my_prop) or spec[i]["c1"].contains(my_prop) {
        totalChance *= chance
    }
)
```

The lower the `totalChance`, the more likely it is the property will change by any modification to the DNA,
 and the less certain we are about the stability of the estimations for this dna.


#### Chance of getting rid of an attribute


```
spec = outputJSON array

totalChance = 1

for each DNA bit "i" [range [0...255]] (
    var chanceP;
    if currentDNA_P[i] == 1 (
        chanceP = spec[i]["c1"][my_prop] || 0     <--- 0 if my_prop key is absent
    ) else (
        chanceP = spec[i]["c0"][my_prop] || 0
    )
    var chanceQ;
    var chanceX;
    if currentDNA_Q[i] == 1 (
        chanceQ = spec[i]["c1dom"][my_prop] || 0
        chanceX = spec[i]["c1dom"][my_prop]
    ) else (
        chanceQ = spec[i]["c0dom"][my_prop] || 0
        chanceX = spec[i]["c0dom"][my_prop]
    )
    var maxCount = spec[i]["max"]
    # Normalize the counts, now they're chances, range [0...1]
    chanceP = chanceP / maxCount
    chanceQ = chanceQ / maxCount
    chanceX = chanceX / maxCount
    
    # This is just the chance for this single bit contributing to the property not being expressed anymore
    # [The chance of the property not showing up after a bitflip] =
    #  [chance of the dominant bit disappearing and not having a recessive expression] + [chance of chance to other dominant bit]
    var chance = ((1 - chanceP) * (1 - chanceX)) + (chanceX - chanceQ)
    
    // If there's no chance for either 0 or 1, recessive (string P),
    //  then this bit does not affect the expression of the property.
    if spec[i]["c0"].contains(my_prop) or spec[i]["c1"].contains(my_prop) {
        // chance is cumulative here
        totalChance += chance
    }
)

Note: the chances of each bit randomly mutating are not strictly independent. And there's crossovers as well, which are even more dependent chances. 
```


#### Chance of getting an attribute *with known `P` and `Q`**

```
spec = outputJSON array

totalChance = 1

for each DNA bit "i" [range [0...255]] (
    var chanceP;
    if currentDNA_P[i] == 1 (
        chanceP = spec[i]["c1"][my_prop] || 0     <--- 0 if my_prop key is absent
    ) else (
        chanceP = spec[i]["c0"][my_prop] || 0
    )
    var chanceQ;
    var chanceX;
    if currentDNA_Q[i] == 1 (
        chanceQ = spec[i]["c1dom"][my_prop] || 0
        chanceX = spec[i]["c1dom"][my_prop]
    ) else (
        chanceQ = spec[i]["c0dom"][my_prop] || 0
        chanceX = spec[i]["c0dom"][my_prop]
    )
    var maxCount = spec[i]["max"]
    # Normalize the counts, now they're chances, range [0...1]
    chanceP = chanceP / maxCount
    chanceQ = chanceQ / maxCount
    chanceX = chanceX / maxCount
    
    # This is just the chance for this single bit contributing to the property being expressed
    # chance of string P getting expressed and having the property + chance of Q being dominant and having it
    var chance = (chanceP * (1 - chanceX)) + chanceQ
    
    // If there's no chance for either 0 or 1, recessive (string P),
    //  then this bit does not affect the expression of the property.
    if spec[i]["c0"].contains(my_prop) or spec[i]["c1"].contains(my_prop) {
        totalChance *= chance
    }
)

```

#### Chance of acquiring an attribute, within the fusion of two chips i.e. semi-unknown P and Q

```

chanceFn (P string, Q string) -> from above, Chance of getting a property [with known `P` and `Q`]

# The parent may transition one of their strings to the child. (simplified, crossovers are involved too)
chance =
   (1/4 * chanceFn(fatherP, motherP))
 + (1/4 * chanceFn(fatherP, motherQ))
 + (1/4 * chanceFn(fatherQ, motherP))
 + (1/4 * chanceFn(fatherQ, motherQ)) + 
```

