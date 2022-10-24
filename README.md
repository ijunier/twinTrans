# twinTrans
Python implementation of the twin transcriptional-loop model using a Gillespie algorithm

---

## Usage

The code has been tested using `Python 3.7`.

### Basic

The most basic usage consists in running:

```
bin/twin.py {results_directory}
```

This will run a simulation with default parameters (lasting $\sim$ 2 min on a 3.1 GHz Intel Core i7). The outcome is composed of two files written out in `{results_directory}`: 

- **param_var.txt**: parameters and variables (and their value) of the simulation 
- **mean_properties.txt**: values of various quantities of interest (written out every a fixed number of transcripts as specified by `-Net` optional argument):
  - *transcripts_nb*: number of transcripts
  - *time*: real time ($s$)
  - *prod_rate*: production rate (= transcripts_nb/time) ($s^{-1}$)
  - *mean_prod_time*: average time separating two successive productions of a transcript (should be equal to 1/prod_rate) ($s$)
  - *mean_bind_time*: average time separating two successive binding events ($s$)
  - *mean_ocf_time*: average time to form the open complex once the RNAP is bound at the promoter ($s$)
  - *mean_esc_time*: average time to escape the promoter once the open complex is formed ($s$)
  - *mean_init_time*: average time between two successive initiations of elongation ($s$)
  - *mean_elong_time*: average elongation time ($s$)

To specificy parameters such as those associated with the promoter ($k_b$, $k_o$, $\sigma_o$ and $k_e$), run:

    bin/twin.py -h

### Promoter-following mode

To follow the topological properties at the promoter, use the option `-promfollow`


    
