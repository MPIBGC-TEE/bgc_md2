#!/bin/bash
# To avoid unintentional removal of subrepos do this in an extra directory where you have clondes bgc_md2 (typically test branch) without recursion to submodules
branch_name="binder_orphan_no_subrepos"

parent=$(git branch --show-current)
# remove old remote and local versions of the branch  if they exist
if git ls-remote --quiet --exit-code --heads origin ${branch_name}  
then 
  echo "found remote branch ${branch_name} and will remove it"  ;  git push -d origin ${branch_name}
else
  echo "remote branch ${branch_name} does not exist."
fi

if git show-ref  --quiet refs/heads/${branch_name}
then 
  echo "erasing local branch ${branch_name}."
  git branch -D ${branch_name}
else
  echo "local branch ${branch_name} does not exist."
fi

# create new branch
git checkout --orphan ${branch_name}
# remove the large directories
rm .gitmodules
rm -rf talks
rm -rf posters
rm -rf manuscripts
git commit -am "automatically created by ${0}  branch ${branch_name} from branch ${parent}. This branch has .gitmodules removed
and has been created without history to provide a minimial source for installation directly via pip from github."
git push --set-upstream origin ${branch_name}

# go back to the original branch
git checkout ${parent}
