# transformer-bio

# How to contribute?

1. Make sure that you have the latest version of the code and that you are standing on the branch "main".

	-`git branch` <- to check where are you standing
	-`git checkout main` <- to move to branch "main"
	-`git pull`

2. Create a branch to write all your code. Make sure that you're on the new branch before writing anything. **The name of your branch should always begin with your name** (i.e. 'alejandro/some-name-for-the-branch-related-to-a-task')

	-`git checkout -b 'your-name/name-of-the-branch'` <- to create and move to the new branch.
	-`git branch` <- to check that you're actually in the right place. If not, got to step 1.

	Note: it's good practice to create a branch for every task you'll be coding. Do not put several tasks in one branch please, it will get messy.

1. Write your code and push it to your branch. Whilst on your branch, write all your code. Once you're done, do the following:

	-`git status` <- to check your changes.
	-`git add what-you-want-to-add`  <- add the specific files you want to add
	-`git commit -m` "some reference message" <- in case we need to go back
	-`git push origin name-of-your-branch`

4. Merge your branch with "main". Now that all your changes are correct, move to "main" and merge.

	-`git checkout main` <- to move to "main"
	-`git pull` <- make sure to have the latest version of the code.
	-`git merge name-of-your-branch` <- this could result in some conflicts regarding code consistency.

	Note: solve code inconsistencies if needed. Your editor should provide a nice UI to do this, or use the terminal if you're hardcore.

5. Create a pull request. Once you do the merging, go to the link that comes out in the terminal (or go to the repo directly on GitHub) and create a pull request. There you can put a description of what it is that your task is doing and assign someone for code review.

6. Finally, send the pull request to the person that's going to do the code review and hope that you didn't break anything.
