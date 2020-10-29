Barnosky_score <- function(x) {

  score <- integer(nrow(x))
  method <- integer(nrow(x))
  Stratigraphy <- integer(nrow(x))
  material <- integer(nrow(x))
  archaeo <- integer(nrow(x))

  for (i in 1:nrow(x)) {
    a <- x[i, ]

    score.me <- NA
    score.as <- NA
    score.ma <- NA
    score.ar <- NA

    # Method

    if (a$Dating_Method == "No info") {
      score.me <- 1
    } else if (str_detect(a$Dating_Method, "AMS")) {
      score.me <- 2
    } else {
      score.me <- 1
    }

    # Stratigraphy
    if (a$Date_type == "Direct") {
      score.as <- 5
    } else if (a$Date_type == "Indirect") {
      #if (a$Stratigraphy == "a" | a$Stratigraphy == "b") {
        #score.as <- 2
      #}

      if (Stratigraphy == "No info" | a$Archaeo == "Associated date") {
        score.as <- 1
      }

      if (a$Stratigraphy == "Same stratum") {
        score.as <- 3
      } else { score.as <- 1 }
    }

    # Material

    if (a$Material_dated == "No info") {
      score.as <- 0
    } else if (a$Material_dated == "Collagen" | a$Material_dated == "Feces" | a$Material_dated == "Dung" | a$Material_dated == "Coprolite" |
        a$Material_dated == "Skin" |a$Material_dated == "Flesh" | a$Material_dated == "Hair" |
        a$Material_dated == "Hydroxyproline" | a$Material_dated == "Keratin") {
      score.ma <- 5
    } else if (a$Material_dated == "Bone apatite") {
      score.ma <- 3
    } else if (a$Material_dated == "Bone" | a$Material_dated == "Mixed bones" |
        a$Material_dated == "Terrestrial carbonate" | a$Material_dated == "Carbonate") {
      score.ma <- 1
    } else if (a$Material_dated == "Charcoal") {
      score.ma <- 6
    } else if (a$Material_dated == "Wood" | a$Material_dated == "alpha-cellulose") {
      score.ma <- 5
    } else if (a$Material_dated == "Peat" | a$Material_dated == "Organic soil" |
        a$Material_dated == "Mud") {
      score.ma <- 3
    } else if (a$Material_dated == "Shell") {
      score.ma <- 2
    } else {
      score.ma <- 0
      }

    # Archaeology
    if(a$Archaeo == "No info") {
      score.ar <-  0
    } else if (a$Archaeo == "assemblage") {
      score.ar <- 4
    } else if (a$Archaeo == "Clear hearth" | a$Archaeo == "Single artifact") {
      score.ar <- 3
    } else if (a$Archaeo == "Clear butchering") {
      score.ar <- 2
    } else if (a$Archaeo == "human remains") {
      score.ar <- 5
    } else {
      score.ar <- 0
    }

    # Assign values
    method[i] <- score.me
    Stratigraphy[i] <- score.as
    material[i] <- score.ma
    archaeo[i] <- score.ar
    score[i] <- score.me + score.as + score.ma + score.ar
  }

  result <- data.frame(AgeID = x$lab.no,
                       Method = method, Stratigraphy = Stratigraphy, Material = material, Archaeo = archaeo,
                       Final.Score = score, ID = x$ID, stringsAsFactors = FALSE)

  # Create column indicating if a given score can be considered reliable or not
  result$Quality <- ifelse(result$Final.Score > 10, "Reliable", "Unreliable")
  result$Quality[which(is.na(result$Quality))] <- "Insufficient information to assess"

  return(result)
}
