# Must be sourced **AFTER** "./R/utils-0_project_settings.R"

# on RSPH: cpu-per-task = 8, but only 7 at a time, mem-per-cpu = 5G
hpc_configs <- EpiModelHPC::swf_configs_rsph(
  partition = "preemptable",
  r_version = "4.3.0",
  mail_user = mail_user
)
#
# # on mox: 20 sim per node, full mem, 24h
# hpc_configs <- EpiModelHPC::swf_configs_hyak(
#   hpc = "mox",
#   partition = "ckpt",
#   r_version = "4.1.2",
#   mail_user = mail_user
# )

# hpc_configs <- EpiModelHPC::swf_configs_hyak(
#   hpc = "klone",
#   partition = "ckpt",
#   r_version = "4.1.1",
#   mail_user = mail_user
# )
