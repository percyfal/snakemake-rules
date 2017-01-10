#! /bin/bash
#
# Copyright (C) 2016 by Per Unneberg
#
# Description: automate build and deployment
# 

# Script directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Variables
RECIPE_DIR=$DIR/../conda.recipe
CONDAMETA=$RECIPE_DIR/meta.yaml
# Current branch
current_branch=$(git rev-parse --abbrev-ref HEAD)

# a simple help
if [ "$1" == "-h" ]; then
    usage="$(basename "$0") [-h] -- program to trigger the deployment of releases and devel pkgs
    where:
        -h     show this help text
        -r     new release tag in the form X.X.X or X.X-(alpha|beta).X
        -d     devel build tag in the form X.X.X[dev]number, e.g. 0.1.0dev1 or X.X-(alpha|beta).X[dev]number, e.g. 0.1-alpha.1dev1
    "
    echo "$usage"
    exit 0
fi

# option management
while getopts r:d: option
do
    case "${option}" in
        r) release=${OPTARG};;
	d) devel=${OPTARG};;
    esac 
done

# release tags
release_short=$( echo $release | sed 's|\([^\^.]*\)\(\^0\)$|\1|g' | sed 's/-alpha/a/g' | sed 's/-beta/b/g' )
release_branch=release_$release

# main
if [[ -z "$devel" && ! -z "$release" ]]; then
    echo "You have triggered the release process"

    # First check that tag doesn't exist
    if git rev-parse $release > /dev/null 2>&1
    then
	echo "Tag $release already exists; aborting release"
	exit 1
    fi
    # Then check that release doesn't exist
    if git rev-parse $release_branch > /dev/null 2>&1
    then
	echo "Release branch $release_branch already exists; aborting release"
	exit 1
    fi
    
    # create a new branch
    git checkout -b $release_branch

    # Add version tag
    git tag -a $release -m "Version $release"

    # CHANGELOG generation
    gitchangelog > $DIR/../ChangeLog
    ret=$?
    if [ ! -e $DIR/../sphinx/source/docs/releases/$release.rst ]; then ret=1; fi;
    if [ ! `grep $release $DIR/../sphinx/source/docs/release_notes.rst` ]; then ret=1; fi;
    if [ $ret -ne 0 ]; then
        echo "Exiting because ChangeLog generation failed."
        echo "Check you actually have a $DIR/../sphinx/source/docs/releases/<tag>/.rst file."
        echo "And you actually updated $DIR/../sphinx/source/docs/release_notes.rst with the latest releases/<tag>."
        exit 1
    fi
    git add $DIR/../ChangeLog
    git commit -m "Updating ChangeLog."

    # conda version number updates
    sed -i "s/release = .*/release = \"$release\" %}/g" $CONDAMETA
    sed -i "s/release_short = .*/release_short = \"$release_short\" %}/g" $CONDAMETA
    git add $CONDAMETA
    echo "Wrote $CONDAMETA"
    git commit -m "Updating conda version to $release ($release_short)."

    # Merge branch into master and push to origin
    git checkout master
    git pull origin
    git merge --no-ff $release_branch -m "Merge branch $release_branch"
    git push origin master
    git branch -d $release_branch

    git tag -a $release -f
    git push origin $release

    git checkout develop
    git merge master
    git push origin develop
    make -f sphinx/Makefile gh-pages

    # Trigger conda build
    conda build conda.recipe
    CONDA_OUTPUT=`conda build conda.recipe/ --output`
    if [ -e $CONDA_OUTPUT ]; then
	anaconda upload $CONDA_OUTPUT
    else
	echo "Failed to upload $CONDA_OUTPUT: $?"
    fi

elif [[ ! -z "$devel" && -z "$release" ]]; then
    echo "You have triggered the devel build process"

    # Check the tag
    if git rev-parse $devel > /dev/null 2>&1
    then
	echo "Tag $devel already exists; aborting devel build"
	exit 1
    fi

    # Push the current branch; should be develop or a feature branch
    if [[ ! $current_branch == feature* ]] && [[ ! $current_branch == develop ]]; then
	echo "Branch '$current_branch' is neither the develop or a feature branch; aborting devel build"
	exit 1
    else
	git push origin $current_branch
    fi

    # tag it locally
    git tag -a $devel -m "New devel[rc] build $devel."

    # and push the tag
    git push origin $devel
    echo "The new devel build was triggered."
    make -f sphinx/Makefile gh-pages-dev

else
    echo "You have to pass a -d tag (dev build) OR -r tag for release."
    echo "Run ./deploy.sh -h to get some more help with the args to pass."
    exit 0
fi
