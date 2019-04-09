arg="[()0-9a-zA-Z+*._> -]+"
type="[()a-zA-Z0-9_ *]+"
if false ;then
    #minus sign must be the last or first inside a []. no need to escape ><
    #reformat calloc(s, sizeof(b)) by mycalloc(s, b)
    sed -i "" -E "s|(calloc\($arg),[ ]*sizeof\(($type)\)\);|my\1,\2\);|g" $@
    #reformat malloc(s*sizeof(b) by mymalloc(s, b)
    sed -i "" -E "s|(malloc\($arg)\*[ ]*sizeof\(($type)\)\);|my\1,\2);|g" $@
    #reformat malloc(sizeof(b)*s by mymalloc(s, b)
    sed -i "" -E "s|(malloc)\(sizeof\(($type)\)[ ]*\*[ ]*($arg)\);|my\1(\3,\2);|g" $@
    #reformat realloc(p, s*size(b)) by myrealloc(p, s, b)
    sed -i "" -E "s|(realloc\($arg),($arg)\*[ ]*sizeof\(($type)\)[ ]*\);|my\1,\2,\3);|g" $@
    #reformat realloc(p, sizeof(b)*s by myrealloc(p, s, b)
    sed -i "" -E "s|(realloc\($arg),[ ]*sizeof\(($type)\)[ ]*\*[ ]*($arg)[ ]*\);|my\1,\3,\2);|g" $@
    #reformat realloc(p, s) by myrealloc(p, s, b)
    sed -i "" -E "s|(realloc\($arg),($arg)\);|my\1,\2,char);|g" $@
    #remove type case
    sed -i -E "s/(\($type\))(mymalloc|myrealloc|mycalloc)/\2/g" $@
fi
fun="[a-zA-Z0-9_>.\-]+"
sed -i -E "s|($fun)\[($arg)\]\[($arg)\]|\1\[(\2)*NN+\3\])|g" $@
