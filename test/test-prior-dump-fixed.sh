#Exit on error
set -e

check_equal() {
    if diff $1 $2; then
        return 0
    else
        return 1
    fi
}

./bcall prior-dump-fixed ../test/data/samples.txt temp.dump ../test/data/regions.bed \
        2>prior-dump-fixed-1.stderr 1>prior-dump-fixed-1.stdout

check_equal prior-dump-fixed-1.stdout ../test/data/prior-dump-fixed-expected1.stdout
check_equal prior-dump-fixed-1.stderr ../test/data/prior-dump-fixed-expected1.stderr
