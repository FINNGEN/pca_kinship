

nosetests . > nose.output 2>&1
RET=$?

if [ $RET -ne 0 -o 1 ] ; then
    echo fail

    #mailx -n -S smtp=smtp.aalto.fi -S from='rkd+noreply@zgib.net' \
    #    -s "verkko build failures" \
    #    rkd@zgib.net <<EOF
    cat <<EOF
The verkko tests have failed.


============
`cat nose.output`
EOF
fi