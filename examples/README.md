# Examples

| script | needs | does |
|---|---|---|
| `toy_end_to_end.py` | nothing but the package | writes a 400 bp chromosome, a four-copy family and a bigWig, runs the whole chain, and checks the result against arithmetic |
| `verify_overlay.py` | your own coordinates, genome and bigWig | rebuilds an overlay and compares a sample of copies against a direct reading of the same file |

```bash
python examples/toy_end_to_end.py /tmp/ma_mapper_example

python examples/verify_overlay.py --coordinates copies.bed \
    --genome hg38.fa --bigwig hg38.phyloP447way.bw --sample 10
```

`toy_end_to_end.py` is the shortest complete use of the API: coordinates →
sequence → alignment → overlay → one value per column. It is also the quickest
way to tell whether an installation works, since it checks its own output.

`verify_overlay.py` is the check worth running on a new dataset before building
anything on it. A copy whose independent reading disagrees means either a bug
or, far more often, an alignment and a coordinate table built from different
copy sets.
