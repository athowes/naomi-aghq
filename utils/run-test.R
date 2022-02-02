id <- orderly::orderly_run("prev-anc-art_sim")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_model0")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_model1")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_model2")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_model3")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_model4")
orderly::orderly_commit(id)

id <- orderly::orderly_run("prev-anc-art_process-results")
orderly::orderly_commit(id)
