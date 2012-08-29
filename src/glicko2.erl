-module(glicko2).

-export([configuration/3,
         read_config/1]).

-export([phi_star/2,
         rate/4, rate/5]).

-ifdef(TEST).
-export([glicko_volatility_test/0]).
-export([glicko_test/0]).
-endif.

-define(EPSILON, 0.000001).
-record(config, { rd, v, tau}).

configuration(IRD, IV, Tau) ->
    #config { rd = IRD,
              v  = IV,
              tau = Tau }.

read_config(#config { rd = RD, v = Sigma, tau = Tau}) ->
    {RD, Sigma, Tau}.

square(V) -> V*V.

%% Step 2
scale(R, RD) ->
    Mu = (R - 1500) / 173.7178,
    Phi = RD / 173.7178,
    {Mu, Phi}.

g(Phi) ->
    1 / (math:sqrt(1 + 3*Phi*Phi / (math:pi() * math:pi()))).


e(Mu, Muj, Phij) ->
    1 / (1 + math:exp(-g(Phij) * (Mu - Muj))).

%% The Opponents list ends up like the following for a given user
%% [{Muj, Phij, g(Phij), E(Mu, Muj, Phij), Sj}]
scale_opponents(Mu, Opponents) ->
    [begin
         {Muj, Phij} = scale(Rj, RDj),
         {Muj, Phij, g(Phij), e(Mu, Muj, Phij), Sj}
     end || {Rj, RDj, Sj} <- Opponents].

%% Step 3
update_rating(Opponents) ->
    1 / (lists:sum([square(GPhij) * EMMP * (1 - EMMP)
                    || {_Muj, _Phij, GPhij, EMMP, _Score} <- Opponents])).

%% Step 4
compute_delta(V, Opponents) ->
    V * lists:sum([GPhij * (Sj - EMMP)
                   || {_Muj, _Phij, GPhij, EMMP, Sj} <- Opponents]).

%% Step 5
vol_f(Phi, V, Delta, A, #config { tau = Tau }) ->
    PHI2 = Phi*Phi,
    fun(X) ->
            EX = math:exp(X),
            D2 = Delta*Delta,
            A2 = (PHI2 + V + EX),
            P2 = (X - A) / (Tau * Tau),
            P1 = (EX * (D2 - PHI2 - V - EX))  / (2*A2*A2),
            P1 - P2
    end.

vol_k(K, F, A, #config { tau = Tau} = Conf) ->
    Const = A - K*math:sqrt(Tau*Tau),
    case F(Const) < 0 of
        true ->
            vol_k(K+1, F, A, Conf);
        false ->
            Const
    end.

i_compute_volatility(Sigma, Phi, V, Delta, #config { tau = Tau } = Conf) ->
    A = math:log(Sigma*Sigma),
    F = vol_f(Phi, V, Delta, A, Conf),
    B = case Delta*Delta > Phi*Phi + V of
            true ->
                math:log(Delta*Delta - Phi*Phi - V);
            false ->
                vol_k(1, F, A, Conf)
        end,
    FA = F(A),
    FB = F(B),
    try
        compute_volatility(A, B, F, FA, FB, 100)
    catch
        throw:{iterations_exceeded, _Vals} ->
            lager:error("Error in vol comp: ~p", [{Sigma, Phi, V, Delta, Tau}]),
            exit(bad_vol_comp)
    end.

sign(X) when X > 0 -> 1;
sign(X) when X < 0 -> -1;
sign(0.0)          ->  0.

compute_volatility(A, B, F, FA, FB, 0) ->
    throw({iterations_exceeded, {A, B, F, FA, FB}});
compute_volatility(A, B, _F, _FA, _FB, _) when abs(B - A) =< ?EPSILON ->
    math:exp(A/2);
compute_volatility(A, B, F, FA, FB, K) ->
    %% C is the midpoint:
    C = (A + B) * 0.5, FC = F(C),
    D = C + (C - A) * (sign(FA - FB) * FC) / math:sqrt(FC*FC - FA*FB),
    FD = F(D),
    case sign(FD) /= sign(FC) of
        true ->
            compute_volatility(C, D, F, FC, FD, K-1);
        false ->
            case sign(FD) /= sign(FA) of
                true ->
                    compute_volatility(A, D, F, FA, FD, K-1);
                false ->
                    true = sign(FD) /= sign(FB),
                    compute_volatility(D, B, F, FD, FB, K-1)
            end
    end.

%% Step 6
phi_star(SigmaP, Phi) ->
    math:sqrt(square(Phi) + square(SigmaP)).

% Step 7
new_rating(PhiStar, Mu, V, Opponents) ->
    PhiP = 1 / math:sqrt(
                  (1 / square(PhiStar))
                 + (1 / V)),
    L = [GPhij * (Sj - EMMP)
         || {_Muj, _Phij, GPhij, EMMP, Sj} <- Opponents],
    MuP  = Mu + square(PhiP) * lists:sum(L),
    {MuP, PhiP}.

% Step 8
unscale(MuP, PhiP) ->
    RP = 173.7178*MuP + 1500,
    RDP = 173.7178*PhiP,
    {RP, RDP}.

rate(R, RD, Sigma, Opponents) ->
    rate(R, RD, Sigma, Opponents, configuration(350, 0.06, 0.5)).

rate(R, RD, Sigma, Opponents, Conf) ->
    {Mu, Phi} = scale(R, RD),
    ScaledOpponents = scale_opponents(Mu, Opponents),
    V = update_rating(ScaledOpponents),
    Delta = compute_delta(V, ScaledOpponents),
    SigmaP = i_compute_volatility(Sigma, Phi, V, Delta, Conf),
    PhiStar = phi_star(SigmaP, Phi),
    {MuP, PhiP} = new_rating(PhiStar, Mu, V, ScaledOpponents),
    {R1, RD1} = unscale(MuP, PhiP),
    {R1, RD1, SigmaP}.


-ifdef(TEST).

data() ->
    Player = {a, 1500, 200},
    Volatility = 0.06,
    Opponents = [{1400, 30,  1},
                 {1550, 100, 0},
                 {1700, 300, 0}],
    {Player, Volatility, Opponents}.

within(X, Y) -> abs(X - Y) < 0.0001.

glicko_test() ->
    {{a, R, RD}, Sigma, Opponents} = data(),
    {Mu, Phi} = scale(R, RD),
    {0.0, 1.1512924985234674} = {Mu, Phi},
    ScaledOpponents = scale_opponents(Mu, Opponents),
    [{-0.5756462492617337,0.1726938747785201,
      0.9954980064506083,0.6394677305521533,1},
     
     {0.28782312463086684,0.5756462492617337,
      0.9531489778689763,0.4318423561076679,0},

     {1.1512924985234674,1.726938747785201,
      0.7242354780877526,0.30284072909521925,0}] = ScaledOpponents,
    V = update_rating(ScaledOpponents),
    1.7789770897239976 = V,
    Delta = compute_delta(V, ScaledOpponents),
    -0.4839332609836549 = Delta,
    SigmaP = i_compute_volatility(Sigma, Phi, V, Delta, configuration(350, 0.06, 0.5)),
    true = within(0.059995984286488495, SigmaP),
    PhiStar = phi_star(SigmaP, Phi),
    true = within(1.1528546895801364, PhiStar),
    {MuP, PhiP} = new_rating(PhiStar, Mu, V, ScaledOpponents),
    true = within(-0.20694096667525494, MuP),
    true = within(0.8721991881307343, PhiP),
    {R1, RD1} = unscale(MuP, PhiP),
    true = within(1464.0506705393013, R1),
    true = within(151.51652412385727, RD1).

glicko_volatility_test() ->
    V = 1.7785,
    Delta = -0.4834,
    Tau = 0.5,
    Sigma = 0.06,
    Phi = 200 / 173.7178,
    SigmaP = i_compute_volatility(Sigma, Phi, V, Delta, configuration(350, 0.06, Tau)),
    SigmaP.

-endif.
