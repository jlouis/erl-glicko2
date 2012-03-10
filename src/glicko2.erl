-module(glicko2).

-export([phi_star/2,
         glicko_test/0,
         rate/4]).

-define(TAU, 0.5). % Good values are between 0.3 and 1.2
-define(EPSILON, 0.000001).

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
vol_f(Phi, V, Delta, A) ->
    PHI2 = Phi*Phi,
    fun(X) ->
            EX = math:exp(X),
            D2 = Delta*Delta,
            A2 = (PHI2 + V + EX),
            P2 = (X - A) / (?TAU * ?TAU),
            P1 = (EX * (D2 - PHI2 - V - EX))  / (2*A2*A2),
            P1 - P2
    end.

vol_k(K, F, A) ->
    Const = A - K*math:sqrt(?TAU*?TAU),
    case F(Const) < 0 of
        true ->
            vol_k(K+1, F, A);
        false ->
            Const
    end.

compute_volatility(Sigma, Phi, V, Delta) ->
    A = math:log(Sigma*Sigma),
    F = vol_f(Phi, V, Delta, A),
    B = case Delta*Delta > Phi*Phi + V of
            true ->
                math:log(Delta*Delta - Phi*Phi - V);
            false ->
                vol_k(1, F, A)
        end,
    FA = F(A),
    FB = F(B),
    compute_volatility(A, B, F, FA, FB).

compute_volatility(A, B, _F, _FA, _FB) when abs(B - A) < ?EPSILON ->
    math:exp(A/2);
compute_volatility(A, B, F, FA, FB) ->
    C = A + (A - B)*FA / (FB - FA),
    FC = F(C),
    {NA, NFA} =
        case FC * FB < 0 of
            true ->
                {B, FB};
            false ->
                {A, FA / 2}
        end,
    compute_volatility(NA, C, F, NFA, FC).
            
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
    {Mu, Phi} = scale(R, RD),
    ScaledOpponents = scale_opponents(Mu, Opponents),
    V = update_rating(ScaledOpponents),
    Delta = compute_delta(V, ScaledOpponents),
    SigmaP = compute_volatility(Sigma, Phi, V, Delta),
    PhiStar = phi_star(SigmaP, Phi),
    {MuP, PhiP} = new_rating(PhiStar, Mu, V, ScaledOpponents),
    {R1, RD1} = unscale(MuP, PhiP),
    {R1, RD1, SigmaP}.

data() ->
    Player = {a, 1500, 200},
    Volatility = 0.06,
    Opponents = [{1400, 30,  1},
                 {1550, 100, 0},
                 {1700, 300, 0}],
    {Player, Volatility, Opponents}.

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
    SigmaP = compute_volatility(Sigma, Phi, V, Delta),
    0.059995984286488495 = SigmaP,
    PhiStar = phi_star(SigmaP, Phi),
    1.1528546895801364 = PhiStar,
    {MuP, PhiP} = new_rating(PhiStar, Mu, V, ScaledOpponents),
    {-0.20694096667525494, 0.8721991881307343} = {MuP, PhiP},
    {R1, RD1} = unscale(MuP, PhiP),
    {1464.0506705393013, 151.51652412385727} = {R1, RD1}.

