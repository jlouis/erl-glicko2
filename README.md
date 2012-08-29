# Glicko 2 Ranking code.

This repository contains code for ranking according to Mark E.
Glickmans "Glicko2" ranking system. See

http://glicko.net/glicko.html

# Using the code

There are are couple of exports, which are the main interface to the
code base. Suppose we have

    Rating = 1500,
    RatingDeviation = 200,
    Sigma = 0.06
    Opponents = [{1400, 30,  1},
                 {1550, 100, 0},
                 {1700, 300, 0}].

Where each opponent is given as a triple `{R, RD, S}` of respectively
the rating of the opponent, the rating deviation and the score. A
score of `1` means we won and a score of `0` means we lost. A score of
`0.5` is a draw, but I am not using that in the code.

Calling

    glicko2:rate(Rating, RatingDeviation, Sigma, Opponents).

Returns a triple `{R1, RD1, Sigma1}` of the players new rating, rating
deviation and sigma/volatility value.

Alternatively, you can generate a configuration `configuration(350,
0.06, 0.5)` which defines, respectively, the initial rating deviation,
the standard volatility base and the value of *tau* in Glicko2. The
*tau* value says something about how much volatility means as a
factor. I suggest tuning these for optimal values.

To use your configuration, call:

    Conf = glicko2:configuration(360, 0.06, 0.5),
    glicko2:rate(R, RD, Sigma, Opponents, Conf).

Another export is `phi_star/2` which is an internal step in the
Glicko2 scheme. This step is useful for tuning, hence it is exported.

# License

MIT

