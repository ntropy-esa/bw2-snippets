# Contribution analysis in brightywa2

There are several ways of performing contribution analysis in matrix-based LCA, and several ways of implementing it.

We can distinguish:

1. process contribution, applied to the whole technosphere matrix
2. stressor contribution, applied to the whole biosphere matrix

>> These two rely on "top_contributor" functions. Implemented in bw2 & AB.

3. process contribution, groupped by existing tags on process activities (e.g. name, or user defined tag)
4. stressor contribution, groupped by existing tags on process activities (e.g. name, compartment, regardless of sub-compartment)

>> These two rely on "top_contributor" functions & sums on dataframe rows, based on some criteria. Implemented in AB, easy to impement in bw2 from the lca results.

5. process contribution, on a given foreground system, by projection of background impacts onto foreground systems, optionally followed by grouping by user defined groups. This is the kind of grouping that can be performed in **SimaPro**. A requirement for this contribution analysis to be mathematically valid: the foreground system must not have internal (infinite) loops. It can be performed with matrix operations or graph-traversal algorithms, equivalently. The matrix algebra is present in Strommans LCA course material, while graph-traversal is implemented in bw2 in a specific case, and extended in bw2-contributionAnalysis.py to multiple impact categories / foreground databases.


>> This grouping is partly available in bw2, relying on graph traversal, but is limited to a single foreground database. 
>> It is not implemented in the activity-browser: 

What is needed to add in the activity-browser?
* In the activity editor tab, add a button `Edit other tags`. The button pops an interface to add user-defined tags, as depicted below:
```
"grouping_A": "Group X in classifcation A"
"grouping_B": "Group Z in classifcation B"
[Add New Tag] [Remove Tag] [Edit Tag]
```
* It can be also convenient to tag an entire foreground database, iteratively, as with the function `conveniently_tag_database_v2`

* Then, additional parameters needs to be defined when performing a calculation:
  - how to group biosphere exchanges: bio2tech = True > projects any biosphere exchange in the foreground onto the tag of its technosphere activity
  - how to deal with untagged activities: parent4other = True > untagged activities will get tag of their nearest parent that had a tag, if any, otherwise set to defaultTag (Other)
  - Note: bio2tech & parent4other set to True reproduces results as in SimaPro
  - select which databases to perform graph traversal on: must not be a set of databases with internal loops, otherwise graph cannot terminate

* Apply function `multi_traverse_tagged_databases`

* Process aggregated and non-aggregated graphs, to plot results by group



1. Any other contribution analysis technique ?




