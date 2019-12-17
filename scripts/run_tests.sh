# !/bin/bash

CUR=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd "$DIR"
cd ..

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

    echo -e "testing with default parameters"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > "tests/test_gh.csv"
    python tests/test_gh.py "tests/test_gh.csv" "tests/test_gh.png"

    A0=10

    echo -e "testing with non zero acceleration"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > "tests/test_gh_acc.csv"
    python tests/test_gh.py "tests/test_gh_acc.csv" "tests/test_gh_acc.png"

    A0=0
    X_NOISE=50
    V_NOISE=20

    echo -e "testing with extra noise"
    bin/test_gh "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" "$GH0" "$GH1" "$GH2" > "tests/test_gh_noisy.csv"
    python tests/test_gh.py "tests/test_gh_noisy.csv" "tests/test_gh_noisy.png"
}

function test_kalman1d() {
    echo -e "testing 1 dimensional kalman filter (tracks position and velocity, independently"

    local X_NOISE=5
    local V_NOISE=5
    local A0=0
    local V0=100

    echo -e "testing with default parameters"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > "tests/test_kalman1d.csv"
    python tests/test_kalman1d.py "tests/test_kalman1d.csv" "tests/test_kalman1d.png"

    A0=5

    echo -e "testing with non zero acceleration"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > "tests/test_kalman1d_acc.csv"
    python tests/test_kalman1d.py "tests/test_kalman1d_acc.csv" "tests/test_kalman1d_acc.png"

    A0=0
    X_NOISE=10

    echo -e "testing with extra noise"
    bin/test_kalman1d "$N_STEPS" "$DT" "$X_NOISE" "$V_NOISE" "$A0" "$V0" > "tests/test_kalman1d_noisy.csv"
    python tests/test_kalman1d.py "tests/test_kalman1d_noisy.csv" "tests/test_kalman1d_noisy.png"
}

function test_kalman() {
    echo -e "testing multidimensional kalman filter (tracks x and y positions. Velocities are tracked but unobserved)"

    local V_X=100
    local V_Y=100
    local X_NOISE=5
    local Y_NOISE=5
    local V_X_NOISE=10
    local V_Y_NOISE=10
    local A_X=0
    local A_Y=0
    local Q_VAR=10

    echo -e "testing with default parameters"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$V_X_NOISE" "$V_Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > "tests/test_kalman.csv"
    python tests/test_kalman.py "tests/test_kalman.csv" "tests/test_kalman.png"

    A_X=50
    A_Y=20
    Q_VAR=50

    echo -e "testing with non zero acceleration"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$V_X_NOISE" "$V_Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > "tests/test_kalman_acc.csv"
    python tests/test_kalman.py "tests/test_kalman_acc.csv" "tests/test_kalman_acc.png"

    A_X=0
    A_Y=0
    V_X_NOISE=100
    V_Y_NOISE=50
    Q_VAR=100

    echo -e "testing with extra process noise"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$V_X_NOISE" "$V_Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > "tests/test_kalman_process_noisy.csv"
    python tests/test_kalman.py "tests/test_kalman_process_noisy.csv" "tests/test_kalman_process_noisy.png"

    V_X_NOISE=10
    V_Y_NOISE=10
    Q_VAR=10
    X_NOISE=75
    Y_NOISE=20

    echo -e "testing with extra sensor noise"
    bin/test_kalman "$N_STEPS" "$DT" "$V_X" "$V_Y" "$X_NOISE" "$Y_NOISE" "$V_X_NOISE" "$V_Y_NOISE" "$A_X" "$A_Y" "$Q_VAR" > "tests/test_kalman_sensor_noisy.csv"
    python tests/test_kalman.py "tests/test_kalman_sensor_noisy.csv" "tests/test_kalman_sensor_noisy.png"
}

echo -e "Starting statistical filters tests"

test_gh
test_kalman1d
test_kalman

echo -e "Done"

cd "$CUR"
