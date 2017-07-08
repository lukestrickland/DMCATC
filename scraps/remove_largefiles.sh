git filter-branch --tree-filter 'rm -rf `cat /d/Software/DMC_ATCPMDC/large_files.txt | cut -d " " -f 2`  
'--prune-empty <BRANCHES>