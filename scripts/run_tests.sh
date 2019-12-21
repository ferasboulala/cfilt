# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"/../
pwd
./runbuild.sh

N_STEPS=100
DT=0.1

function test_gh() {
    echo -e "testing gh filter (tracks position, velocity and acceleration)"

    local X_NOISE=20
    local V_NOISE=10
    local A0=0
    local V0=100
    local GH0=0.5
    local GH1=0.5
    local GH2=0.5

    local TDIR="tests/gh/"

    echo -e "testing with default parameters"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > ""$TDIR"test_gh.csv"
    python "$TDIR"test_gh.py ""$TDIR"test_gh.csv" ""$TDIR"test_gh.png"

    A0=10

    echo -e "testing with non zero acceleration"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > ""$TDIR"test_gh_acc.csv"
    python "$TDIR"test_gh.py ""$TDIR"test_gh_acc.csv" ""$TDIR"test_gh_acc.png"

    A0=0
    X_NOISE=50
    V_NOISE=20

    echo -e "testing with extra noise"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > ""$TDIR"test_gh_noisy.csv"
    python "$TDIR"test_gh.py ""$TDIR"test_gh_noisy.csv" ""$TDIR"test_gh_noisy.png"
}

function test_kalman1d() {
    echo -e "testing 1 dimensional kalman filter (tracks position and velocity, independently"

    local X_NOISE=5
    local V_NOISE=5
    local A0=0
    local V0=100

    local TDIR="tests/kalman1d/"

    echo -e "testing with default parameters"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > ""$TDIR"test_kalman1d.csv"
    python "$TDIR"test_kalman1d.py ""$TDIR"test_kalman1d.csv" ""$TDIR"test_kalman1d.png"

    A0=5

    echo -e "testing with non zero acceleration"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > ""$TDIR"test_kalman1d_acc.csv"
    python "$TDIR"test_kalman1d.py ""$TDIR"test_kalman1d_acc.csv" ""$TDIR"test_kalman1d_acc.png"

    A0=0
    X_NOISE=10

    echo -e "testing with extra noise"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > ""$TDIR"test_kalman1d_noisy.csv"
    python "$TDIR"test_kalman1d.py ""$TDIR"test_kalman1d_noisy.csv" ""$TDIR"test_kalman1d_noisy.png"
}

function test_kalman() {
    echo -e "testing multidimensional kalman filter (tracks x and y positions. Velocities are tracked but unobserved)"

    local V_X=10
    local V_Y=10
    local X_NOISE=1
    local Y_NOISE=1
    local A_X=0
    local A_Y=0
    local Q_VAR=1

    local TDIR="tests/kalman/"

    echo -e "testing with default parameters"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > ""$TDIR"test_kalman.csv"
    python "$TDIR"test_kalman.py ""$TDIR"test_kalman.csv" ""$TDIR"test_kalman.png"

    A_X=1
    A_Y=2
    Q_VAR=5

    echo -e "testing with non zero acceleration"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > ""$TDIR"test_kalman_acc.csv"
    python "$TDIR"test_kalman.py ""$TDIR"test_kalman_acc.csv" ""$TDIR"test_kalman_acc.png"

    A_X=0
    A_Y=0
    Q_VAR=1
    X_NOISE=5
    Y_NOISE=10

    echo -e "testing with extra sensor noise"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > ""$TDIR"test_kalman_sensor_noisy.csv"
    python "$TDIR"test_kalman.py ""$TDIR"test_kalman_sensor_noisy.csv" ""$TDIR"test_kalman_sensor_noisy.png"
}

echo -e "Starting statistical filters tests"

test_gh
test_kalman1d
test_kalman

echo -e "Done"

cd "$CUR"
