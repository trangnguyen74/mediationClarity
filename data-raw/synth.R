library(tidyverse)


### GET DATA ###################################################################

dat <- foreign::read.spss(here::here("data-raw", "06070809_PAS_data.sav"),
                          use.value.labels = FALSE,
                          to.data.frame = TRUE)

### Copy names of measurement items of mediators at four time points from codebook file: Edit1_selectie_paper Trang.sps, including

# self-control (adolescent)
tmp <- c(paste0(c("a", "b", "e", "f", "g", "h", "i", "k", "l"),
                "x"),
         "c", "d", "j", "m")
controlitems0 <- paste0("xv33", tmp)
controlitems1 <- paste0("yv33", tmp)
controlitems2 <- paste0("zv33", tmp)
controlitems3 <- paste0("wv33", tmp)

# attitudes (adolescent)
tmp <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
attitudeitems0 <- paste0("xv29", tmp)
attitudeitems1 <- paste0("yv29", tmp)
attitudeitems2 <- paste0("zv29", tmp)
attitudeitems3 <- paste0("wv29", tmp)

# adolescent-reported (alcohol related) rules
tmp <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
rulesitems0 <- paste0("xv30", tmp)
rulesitems1 <- paste0("yv30", tmp)
rulesitems2 <- paste0("zv30", tmp)
rulesitems3 <- paste0("wv30", tmp)

# parent attitudes
tmp <- c("a", "b", "c", "d", "e", "f", "g", "h")
paattitudeitems0 <- paste0("av23", tmp)
paattitudeitems1 <- paste0("bv23", tmp)
paattitudeitems2 <- paste0("cv23", tmp)
paattitudeitems3 <- paste0("dv23", tmp)

# parent-reported rules
tmp <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
parulesitems0 <- paste0("av25", tmp)
parulesitems1 <- paste0("bv25", tmp)
parulesitems2 <- paste0("cv25", tmp)
parulesitems3 <- paste0("dv25", tmp)

rm(tmp)



### CLEAN AND CODE VARIABLES ###################################################

rename_each <- function(df, fun){ setNames(df, fun(names(df))) }

dat <- dat %>%
    rename_each(tolower) %>%

    select(-c(school, class)) %>%

    # grab 2 conditions: combined intervention and usual care
    filter(conditie > 2) %>%
    mutate(treat = 1 * (conditie==3),
           treat = factor(treat)) %>%
    select(-conditie) %>%

    rename(sex = "gender",
           edu = "leveleduc",
           age = "aget1") %>%

    # recode one 10-year-olde to age 11 and eighteen 14-year-olds to age 13
    mutate(age = ifelse(age==10, 11,
                        ifelse(age==14, 13,
                               age))) %>%
    mutate(age = factor(age)) %>%

    # recode religion 88 to 5 (other religion)
    mutate(religion = ifelse(religion==88, 5, religion)) %>%

    # parent education as higher level of two parents
    rename(paedu = "leveleduc1") %>%
    mutate(paedu = ifelse(is.na(paedu),
                          leveleduc2,
                          ifelse((!is.na(leveleduc2) & leveleduc2 > paedu),
                                 leveleduc2,
                                 paedu))) %>%
    select(-leveleduc2) %>%

    # weekly drinking at 4 time points (dichotomize from frequency b/c largely zeros)
    rename(wkdrink0 = "qftot",
           wkdrink1 = "qftot_y",
           wkdrink2 = "qftot_z",
           wkdrink3 = "qftot_w") %>%
    mutate(wkdrink0 = 1 * (wkdrink0 > 0),
           wkdrink1 = 1 * (wkdrink1 > 0),
           wkdrink2 = 1 * (wkdrink2 > 0),
           wkdrink3 = 1 * (wkdrink3 > 0)) %>%

    # self-control
    rename(zv33ix = "zv33i") %>%  # (fix item that not yet reverse coded)
    mutate(zv33ix = ifelse(zv33ix==99,
                           zv33ix,
                           6 - zv33ix)) %>%
    mutate_at(c(controlitems0, controlitems1, controlitems2, controlitems3),
              list(~na_if(., 99))) %>%
    mutate(sfc0 = rowMeans(select(., !!controlitems0), na.rm = TRUE),
           sfc1 = rowMeans(select(., !!controlitems1), na.rm = TRUE),
           sfc2 = rowMeans(select(., !!controlitems2), na.rm = TRUE),
           sfc3 = rowMeans(select(., !!controlitems3), na.rm = TRUE)) %>%

    # attitudes (adolescent)
    mutate_at(c(attitudeitems0, attitudeitems1, attitudeitems2, attitudeitems3),
              list(~na_if(., 99))) %>%
    mutate(att0 = rowMeans(select(., !!attitudeitems0), na.rm = TRUE),
           att1 = rowMeans(select(., !!attitudeitems1), na.rm = TRUE),
           att2 = rowMeans(select(., !!attitudeitems2), na.rm = TRUE),
           att3 = rowMeans(select(., !!attitudeitems3), na.rm = TRUE)) %>%

    # adolescent-reported rules
    mutate_at(c(rulesitems0, rulesitems1, rulesitems2, rulesitems3),
              list(~na_if(., 99))) %>%
    mutate(rul0 = rowMeans(select(., !!rulesitems0), na.rm = TRUE),
           rul1 = rowMeans(select(., !!rulesitems1), na.rm = TRUE),
           rul2 = rowMeans(select(., !!rulesitems2), na.rm = TRUE),
           rul3 = rowMeans(select(., !!rulesitems3), na.rm = TRUE)) %>%

    # code non-reponse to wkdrink0 to a separate category (-1)
    mutate(wkdrink0 = ifelse(is.na(wkdrink0), - 1, wkdrink0)) %>%

    # code missing wkdrink1 to -1 if any of rul1, att1, sfc1 is observed
    # (-1 means refuse to answer, not absent)
    # and same with wkdrink2
    mutate(wkdrink1 = ifelse(!is.na(wkdrink1), wkdrink1,
                             ifelse(is.na(rul1) & is.na(att1) & is.na(sfc1), NA,
                                    -1)),
           wkdrink2 = ifelse(!is.na(wkdrink2), wkdrink2,
                             ifelse(is.na(rul2) & is.na(att2) & is.na(sfc2), NA,
                                    -1))) %>%

    # parent attitudes
    mutate_at(c(paattitudeitems0, paattitudeitems1, paattitudeitems2,
                paattitudeitems3),
              list(~na_if(., 99))) %>%
    mutate(paatt0 = rowMeans(select(., !!paattitudeitems0), na.rm = TRUE),
           paatt1 = rowMeans(select(., !!paattitudeitems1), na.rm = TRUE),
           paatt2 = rowMeans(select(., !!paattitudeitems2), na.rm = TRUE),
           paatt3 = rowMeans(select(., !!paattitudeitems3), na.rm = TRUE)) %>%

    # parent-reported rules
    mutate_at(c(parulesitems0, parulesitems1, parulesitems2, parulesitems3),
              list(~na_if(., 99))) %>%
    mutate(parul0 = rowMeans(select(., !!parulesitems0), na.rm = TRUE),
           parul1 = rowMeans(select(., !!parulesitems1), na.rm = TRUE),
           parul2 = rowMeans(select(., !!parulesitems2), na.rm = TRUE),
           parul3 = rowMeans(select(., !!parulesitems3), na.rm = TRUE)) %>%

    # code variables as factor variables
    mutate(sex = factor(sex),
           edu = factor(edu),
           religion = factor(religion),
           paedu = factor(paedu),
           wkdrink0 = factor(wkdrink0),
           wkdrink1 = factor(wkdrink1),
           wkdrink2 = factor(wkdrink2),
           wkdrink3 = factor(wkdrink3)) %>%

    select(code, treat, sex, edu, religion, age, paedu,
           wkdrink0, wkdrink1, wkdrink2, wkdrink3,
           sfc0, sfc1, sfc2, sfc3,
           att0, att1, att2, att3,
           rul0, rul1, rul2, rul3,
           paatt0, paatt1, paatt2, paatt3,
           parul0, parul1, parul2, parul3) %>%

    # leave out parent variables
    select(-starts_with("pa"))




### SINGLE IMPUTATION OF MISSING DATA ##########################################

dat <- dat %>% select(-code)

dat0 <- dat %>% filter(treat==0)
dat1 <- dat %>% filter(treat==1)

# view missing patterns
tmp <- mice::md.pattern(dat0, rotate.names = TRUE)
tmp <- mice::md.pattern(dat1, rotate.names = TRUE)
rm(tmp)


if (exists(".Random.seed", .GlobalEnv)) { oldRNG <- .GlobalEnv$.Random.seed
} else                                  { oldRNG <- NULL
}

imp0 <- mice::mice(dat0, method = "cart", seed = 33, m = 1)
imp1 <- mice::mice(dat1, method = "cart", seed = 77, m = 1)

if (!is.null(oldRNG)) { .GlobalEnv$.Random.seed <- oldRNG
} else                { rm(".Random.seed", envir = .GlobalEnv)
}
rm(oldRNG)


idat <- bind_rows(mice::complete(imp0, 1),
                  mice::complete(imp1, 1))



### CREATE SYNTHETIC DATASET ###################################################

idat <- idat %>%
    select(sex, age, religion, edu,
           wkdrink0, att0, sfc0, rul0,
           treat,
           att1, rul1, sfc1, wkdrink1,
           att2, rul2, sfc2, wkdrink2,
           att3, rul3, sfc3, wkdrink3)


if (exists(".Random.seed", .GlobalEnv)) { oldRNG <- .GlobalEnv$.Random.seed
} else                                  { oldRNG <- NULL
}

synth <- synthpop::syn.strata(idat, strata = "treat", seed = 55)
synth <- synth$syn

if (!is.null(oldRNG)) { .GlobalEnv$.Random.seed <- oldRNG
} else                { rm(".Random.seed", envir = .GlobalEnv)
}
rm(oldRNG)




### ADDITIONAL CODING ##########################################################

synth <- synth %>%
    select(-c(wkdrink1, wkdrink2, rul2, att2, sfc2, rul3, att3, sfc3))


# dichotomize skewed variables
synth %>%
    select(c(sfc0, sfc1, rul0, rul1, att0, att1, treat)) %>%
    gather(-treat, key = "key", value = "value") %>%
    ggplot(aes(x = value, group = treat, color = treat)) +
    facet_wrap(. ~ key, ncol = 2) +
    theme_bw() +
    geom_density()

synth <- synth %>%
    mutate(rul0 = factor(1*(rul0>=4.5)),
           rul1 = factor(1*(rul1>=4.5)),
           att0 = factor(1*(att0>=4.5)),
           att1 = factor(1*(att1>=4.5)))

# recode mediators and outcome as numeric
synth <- synth %>%
    mutate(treat = as.integer(as.character(treat)),
           att = as.integer(as.character(att1)),
           rul = as.integer(as.character(rul1)),
           sfc = sfc1,
           drink = as.integer(as.character(wkdrink3)),
           drink0 = wkdrink0) %>%
    mutate(drink0 = factor(drink0, labels = c("NR", 0, 1))) %>%
    select(treat,
           sex, age, edu, religion,
           att0, rul0, sfc0, drink0,
           att, rul, sfc, drink)

usethis::use_data(synth, overwrite = TRUE)
