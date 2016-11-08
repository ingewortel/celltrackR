
load('ref_tracks.RData')

tracks.from.csv <- read.tracks.csv('ref_tracks.csv', sep = ',')
tracks.from.csv.blank.sep <- read.tracks.csv('ref_tracks_blank_line_sep.csv', sep = ',',
                                             track.sep.blankline = T,
                                             time.column = 1, pos.columns = c(2:4))
tracks.from.df <- as.tracks.data.frame(read.csv('ref_tracks.csv'))

test_that("Tracks are loaded correctly", {
  expect_equivalent(tracks.from.csv, ref)
  expect_equivalent(tracks.from.df, ref)
  expect_equivalent(tracks.from.csv.blank.sep, ref)
} )

test_that("Tracks have correct structure", {
  expect_equal(class(tracks.from.df[[1]]), "matrix")
})
