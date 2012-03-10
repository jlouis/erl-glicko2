-module(glicko2_test).

-include_lib("eqc/include/eqc.hrl").

-compile(export_all).

%% Generators
g_r() ->
    ?LET(I, choose(300000, 2900000),
         I / 1000).

g_rd() ->
    ?LET(I, choose(5000, 350000),
         I / 1000).

g_sigma() ->
    ?LET(I, choose(599, 800),
         I / 10000).

g_player() ->
    {g_r(), g_rd(), g_sigma()}.

g_opponent() ->
    {g_r(), g_rd(), elements([0, 0.5, 1])}.

g_opponents() ->
    non_empty(list(g_opponent())).

prop_glicko_terminates() ->
    ?FORALL({{R, RD, Sigma}, Opponents},
            {g_player(), g_opponents()},
            ?TIMEOUT(
               3000,
               begin
                   {_R1, _RD1, _Sigma1} = glicko2:rate(R, RD, Sigma, Opponents),
                   true
               end)).
