# Contribution analysis in brightywa2

There are many ways of performing contribution analysis in matrix-based LCA, relying on different mathematical techniques, sometimes leading to several equivalent implementations. 

### Main mathematical tools: 
- grouping / summation under conditios (aka masks)
- background to foreground projections
- supply-chain graph traversals (equivalent to projection)
- geometric series expansions


### Way of performing contribution analasys:

1. process contribution, applied to the whole technosphere matrix: the direct impacts of each activity are compared
2. stressor contribution, applied to the whole biosphere matrix: the impact from the total emissions for each stressors are compared. 
   
> Implemented in bw2 & AB. Relies on: sorting the dataframes of technosphere / biosphere results. 

3. process contribution, aggregated by existing tags on process activities (e.g. name, or user defined tag)
4. stressor contribution, aggregated by existing tags on process activities (e.g. name, compartment, regardless of sub-compartment)

> Implemented in AB, easy to impement in bw2 from the lca results. Relies on: summing together some rows of the dataframes of technosphere / biosphere results, based on some conditions. The conditions are usually equality conditions within the meta-data. 

5. process contribution, on a given foreground system, by projection of background impacts onto foreground systems, optionally followed by grouping by user defined groups. This is the kind of grouping that can be performed in **SimaPro**. A requirement for this contribution analysis to be mathematically valid: the foreground system must not have internal (infinite) loops. It can be performed with matrix operations or graph-traversal algorithms, equivalently. 

> This grouping is partly available in bw2, relying on graph traversal, but is limited to a single foreground database. It is not implemented in the activity-browser. The matrix algebra is present in Strommans LCA course material (and implemented in Lisa's Master thesis notebooks), while graph-traversal for multiple impact categories & several foreground databases, is available in bw2-contributionAnalysis.py.

6. Other contribution analysis techniques: impact by tiers, from geometric series expansions; bull-eye visualisation (lcopt) based on multi-level graph traversal; tree visualisations ... ?



### What is needed to add graph-traversal techniques in the activity-browser?

* In the activity editor tab, add a button `Edit other tags`. The button pops an interface to add user-defined tags

* It can be also convenient to tag an entire foreground database, iteratively, as with the function `conveniently_tag_database_v2`

* Then, additional parameters needs to be defined when performing a calculation:
  - how to group biosphere exchanges: bio2tech = True > projects any biosphere exchange in the foreground onto the tag of its technosphere activity
  - how to deal with untagged activities: parent4other = True > untagged activities will get tag of their nearest parent that had a tag, if any, otherwise set to defaultTag (Other)
  - Note: bio2tech & parent4other set to True reproduces results as in SimaPro
  - select which databases to perform graph traversal on: must not be a set of databases with internal loops, otherwise graph cannot terminate!

* Apply function `multi_traverse_tagged_databases`

* Process aggregated and non-aggregated graphs, to plot results









