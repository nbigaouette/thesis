#!/bin/bash

# ScribTeX options
scribtex_username=${scribtex_username-"nbigaouette"}
scribtex_project=${scribtex_project-"thesis"}

branch_to_push="master"
#branch_to_push="revision"

scribtex_branch="scribtex"
scribtex_remote="remotescribtex"

function die()
{
    if [[ "${@}" != "" ]]; then
        echo -e "scribtex_git.sh FAILED: ${@}"
    fi
    git checkout ${branch_to_push}
    exit
}

function usage()
{
    echo "Usage:"
    echo "  ./scribtex_git.sh <action>"
    echo ""
    echo "Possible values for \"action\" can be:"
    echo "  setup     Setup local repository for usage with ScribTeX servers. Do this once."
    echo "  cleanup   Cleanup ScribTeX branch and remote."
    echo "  push      Send modifications to ScribTeX servers."
    echo "  pull      Get modifications from ScribTeX servers."
    echo ""
    echo "Make sure the 'scribtex_username' and 'scribtex_project' variables"
    echo "are have the right value for this project!"
    echo "Current values are:"
    echo "  scribtex_username: ${scribtex_username}"
    echo "  scribtex_project: ${scribtex_project}"
}

function setup_remote()
{
    git checkout ${branch_to_push}                              || die "git checkout ${branch_to_push}"
    # Create a remote called "${scribtex_remote}" pointing to ScribTeX
    git remote add ${scribtex_remote} git@git.scribtex.com:${scribtex_username}/${scribtex_project}.git \
                                                                || die "Adding remote \"${scribtex_remote}\""
    # Create a local branch called "${scribtex_branch}"
    #git branch ${scribtex_branch}                               || die "Creating ${scribtex_branch} branch"
    git branch -u ${scribtex_remote}/${scribtex_branch}         || die "Creating ${scribtex_branch} branch"
    # Fetch data from scribtex' servers
    #git fetch ${scribtex_remote}                                || die "Fetching updated data"
    # Set up the local branch "${scribtex_branch}" to track "${scribtex_remote}" remote's "master" branch so a "git pull" will work.
    #git branch --set-upstream ${scribtex_branch} ${scribtex_remote}/master || die "Tracking upstream"
}

function cleanup()
{
    git branch -d ${scribtex_branch}
    git remote rm ${scribtex_remote}
}

function scribtex_push()
{
    # Detect if any stashes exist.
    nb_stash_before=`git stash list | wc -l`
    git stash                                                   || die "Stashing changes."
    nb_stash_after=`git stash list | wc -l`

    git checkout ${scribtex_branch}                             || die "Checking-out \"${scribtex_branch}\" branch."
    git rebase ${branch_to_push}                                || die "Rebasing branch \"${scribtex_branch}\" on updated \"${branch_to_push}\"."
    git push ${scribtex_remote} ${scribtex_branch}:refs/heads/master || die "Pushing \"${scribtex_branch}\" to ScribTeX site."
    git checkout ${branch_to_push}                              || die "Checking-out \"${branch_to_push}\" branch."

    # If number of stashes changed after calling "git stash", then some changes were
    # stashed. In that case, pop them. If no changes were stashed, don't pop any
    # stash since it could be a user's stash.
    if [[ "${nb_stash_before}" != "${nb_stash_after}" ]]; then
        git stash pop                                           || die "Poping stashed changes."
    fi
}

function scribtex_pull()
{
    git checkout ${scribtex_branch}                             || die "Checking-out \"${scribtex_branch}\" branch."
    git pull                                                    || die "Pulling \"${scribtex_branch}\" branch."
    git checkout ${branch_to_push}                              || die "Checking-out \"${branch_to_push}\" branch."
    git merge ${scribtex_branch}                                || die "Merging \"${scribtex_branch}\" branch on \"${branch_to_push}\" branch."
}

function main()
{
    if [[ ${#@} -ne 1 ]]; then
        usage
        die
    fi

    action="${1}"

    if [[ "${action}" == "setup" ]]; then
        setup_remote
    elif [[ "${action}" == "cleanup" ]]; then
        cleanup
    elif [[ "${action}" == "push" ]]; then
        scribtex_push
    elif [[ "${action}" == "pull" ]]; then
        scribtex_pull
    else
        usage
        die "Unknown action! (${action})"
    fi
}

main ${@}
