PROJECT = glicko2

DEPS = lager
dep_lager = https://github.com/basho/lager.git 2.0.3

ERLC_OPTS = +debug_info '+{parse_transform, lager_transform}'

include erlang.mk
