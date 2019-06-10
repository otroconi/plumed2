#!/bin/bash

extension=no

if [ "$1" = --extension ] ; then
  extension=yes
fi

cd ../
for dir in *
do
  
  if test -f "$dir/module.type"
  then
    if [ "$extension" = no ] ; then
        case "$(cat "$dir/module.type")" in
        (always) echo $dir ;;
        (default-on) test -f $dir.off || echo $dir ;;
        (default-off) test -f $dir.on && echo $dir ;;
        esac
    else
        case "$(cat "$dir/module.type")" in
        (default-on) test -f $dir.off || echo $dir ;;
        (default-off) test -f $dir.on && echo $dir ;;
        esac
    fi
  fi
done
