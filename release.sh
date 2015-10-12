#! /bin/bash

# use this script to make a new release
# just type ./release.sh
# and follow instructions

confirm () {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case $response in
        [yY][eE][sS]|[yY]) 
            true
            ;;
        *)
            false
            ;;
    esac
}

update_changelog() {
  echo $1 $2 $3
  awk -v version="$2" -v shortversion="$3" -v date="$4" '{
    if(version == shortversion ".0" && $1=="Version" && $2==shortversion){
      print "Version " shortversion " (" date ")"
      done=1
    } else if($1=="Version" && $2==version){
      print "Version " version " (" date ")"
      done=1
    } else if(!done && $1=="Unreleased" && $2=="changes"){
      print "Version " version " (" date ")"
    } else print
  }' $1 > $1.tmp
  mv $1.tmp $1
}

VALIDATED=
if [ "$1" = --validated ] ; then
  VALIDATED=yes
elif test -n "$1" ; then
  echo ERROR
  echo "cannot understand $1 option"
  exit 1
fi


echo "*** Only use this script if you want to release a new PLUMED version. ***"
echo "*** Follow instructions below, and use Control-C to exit. ***"

ls src README 1>/dev/null 2>/dev/null ||  {
  echo "Launch from root directory"
  exit 1
}

# Guess version number:
branch="$(git status | grep "On branch" | awk '{print $3}')"
for((i=0;;i++))
do
  guess="${branch#v}.$i"
  if ! git tag | grep -q "^v${guess//./\\.}$" ; then
    echo "Guessed next version: ${guess}"
    break;
  fi
done

case "$branch" in
(v2.?*)
  read -r -p "Type the version number (default: ${guess}): " version
  case "$version" in
  ("")
    version="$guess"
  esac
  ;;
(*)
  read -r -p "Type the version number (e.g. 2.1.3 or 2.2b): " version
  ;;
esac


case "$version" in
(2.?*)
  version="${version#2.}"
  ;;
(*)
  echo "ERROR"
  echo "Use a version in form 2.*"
  exit 1
  ;;
esac

shortversion=$(echo "$version" | sed  's/^\([0-9][0-9]*\).*/\1/' )

version=2.$version
shortversion=2.$shortversion

echo "major version v$shortversion"
echo "minor version v$version"

if test $version == $shortversion ; then
  echo "ERROR"
  echo "please add a patch identified so that tag and branch have different names"
  exit 1
fi

if test "$(git tag -l | awk '{if($1=="'v$version'")print "yes"}')" == yes ; then
  echo "tag v$version already exists"
  exit 1
fi

echo "checking out branch v$shortversion"
git checkout v$shortversion

set -e
if ! test "$VALIDATED" ; then
  update_changelog CHANGES/v$shortversion.txt $version $shortversion "coming soon"
  echo 
  msg="Travis tests for v$version

[makedoc]"
  echo "Now I will add an empty commit and push the result to origin"
  echo "I will use the following commands:"
  echo "***"
  echo "git add CHANGES/v$shortversion.txt"
  echo "git commit --allow-empty -m \"$msg\""
  echo "git push origin v$shortversion"
  echo "***"
  confirm || exit
  git add CHANGES/v$shortversion.txt
  git commit --allow-empty -m "$msg"
  git push origin v$shortversion
  echo
  echo "Now you should go at this link:"
  echo "  http://travis-ci.org/plumed/plumed2/builds"
  echo "and wait for travis to finish the tests"
  echo "In case of failures, fix and repeat the procedure"
  echo "Also check the online manual, that should be here:"
  echo "  http://plumed.github.io/doc-v$shortversion"
  echo "In case of success, relaunch this script as \"./release.sh --validated\""
else
  update_changelog CHANGES/v$shortversion.txt $version $shortversion "$(date '+%b %e, %Y' | sed 's/  / /g')"
  {
    grep  \# VERSION 
    echo $version
  } > VERSION.tmp
  mv VERSION.tmp VERSION
  echo "Here's the new VERSION file:"
  echo "***"
  cat VERSION
  echo "***"
  msg="Release v$version

[makedoc]"
  echo "Now I will add it, prepare a release commit, add a tag named v$version"
  echo "push it to origin and create a tgz file"
  echo "I will use the following commands:"
  echo "***"
  echo "git add VERSION"
  echo "git commit --allow-empty -m \"$msg\""
  echo "git tag v$version"
  echo "git push origin v$version"
  echo "git archive -o plumed-$version.tgz --prefix plumed-$version/ v$version"
  echo "***"
  confirm || exit
  git add VERSION
  git commit --allow-empty -m "$msg"
  git tag v$version
  git push origin v$shortversion v$version
  git archive -o plumed-$version.tgz --prefix plumed-$version/ v$version
  echo
  echo "Done!"
  echo
  echo "Please upload the file plumed-$version.tgz on the download directory"
  echo "Remember to notify the mailing list"
fi






