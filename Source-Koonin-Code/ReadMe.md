# ReadMe File

## Code Source

Downloaded from http://bbs.sciencenet.cn/thread-266127-1-1.html

You can change the name of the files by the following command: (from http://stackoverflow.com/questions/20253584/linux-rename-files-to-uppercase)
```bash
$ for f in * ; do mv -- "$f" "$(tr [:lower:] [:upper:] <<< "$f")" ; done
```

## Files

Follow the description of page 232 from the book.

1. `WarmUp` : Text codes. The example codes from the book.

2. `ComUtil` : Common Utility FORTRAN codes.

3. `Example` and `Project` : Physics codes include all examples and projects

4. `Data` : Data extension for example 5.

5. `ComBlk` : Variable Declaration.
