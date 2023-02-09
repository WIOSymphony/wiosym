# Example of ISO 19115/19139 metadata generated with geometa
# This example is valid according to INSPIRE common metadata requirements

generate_metadata<-function(aRasterFile, pred_name){
  use_packages(c("geometa","lubridate","uuid","stringr"))
  aRaster <- raster(aRasterFile)
  uuid <- UUIDgenerate()
  
  metadata_id <- paste("SGU", uuid, sep = "_")
  
  c_yr <- year(Sys.time())
  c_month <- month(Sys.time())
  c_day <- day(Sys.time())
  
  md <- ISOMetadata$new()
  md$setFileIdentifier(metadata_id)
  # md$setParentIdentifier("my-parent-metadata-identifier")
  md$setCharacterSet(product_character_set)
  md$setLanguage(product_language_code)
  md$setDateStamp(Sys.time())
  # md$setMetadataStandardName("ISO 19115:2003/19139")
  md$setMetadataStandardName(product_metadata_standard_name)
  
  md$setMetadataStandardVersion(product_metadata_standard_version)
  md$setDataSetURI(paste0(product_base_url,product_sub_url, product_title))
  
  # Add contact ----------
  
  rp <- ISOResponsibleParty$new()
  rp$setIndividualName("")
  rp$setOrganisationName(product_org)
  rp$setPositionName(product_contact)
  rp$setRole("")
  contact <- ISOContact$new()
  phone <- ISOTelephone$new()
  phone$setVoice(product_phone)
  # phone$setFacsimile("myfacsimile")
  contact$setPhone(phone)
  address <- ISOAddress$new()
  address$setDeliveryPoint(product_postal1)
  address$setCity(product_city)
  address$setPostalCode(product_postal2)
  address$setCountry(product_country)
  address$setEmail(product_email)
  contact$setAddress(address)
  res <- ISOOnlineResource$new()
  res$setLinkage(product_base_url)
  # res$setName("someresourcename")
  contact$setOnlineResource(res)
  rp$setContactInfo(contact)
  md$addContact(rp)
  
  
  # VectorSpatialRepresentation ---------- 
  if (product_data_type == "vector") {
    vsr <- ISOVectorSpatialRepresentation$new()
    vsr$setTopologyLevel("geometryOnly")
    geomObject <- ISOGeometricObjects$new()
    geomObject$setGeometricObjectType("surface")
    geomObject$setGeometricObjectCount(5L)
    vsr$setGeometricObjects(geomObject)
    md$addSpatialRepresentationInfo(vsr)
  }
  
  if (product_data_type == "grid") {
    vsr <- ISOGridSpatialRepresentation$new()
    vsr$setNumberOfDimensions(2)
    vsr$setTransformationParameterAvailability(TRUE)
    
    dim1 <- ISODimension$new()
    dim1$setName("row")
    dim1$setSize(dim(aRaster)[1])
    dim1$setResolution(ISOMeasure$new(value = xres(aRaster), uom = "m"))
    vsr$addDimension(dim1)
    
    dim2 <- ISODimension$new()
    dim2$setName("column")
    dim2$setSize(dim(aRaster)[2])
    dim2$setResolution(ISOMeasure$new(value = yres(aRaster), uom = "m"))
    vsr$addDimension(dim2)
    vsr$setCellGeometry("area")
    # xml <- vsr$encode()
    md$addSpatialRepresentationInfo(vsr)
  }
  # ReferenceSystem ----------
  rs <- ISOReferenceSystem$new()
  rsId <- ISOReferenceIdentifier$new(code = ref_syst, codeSpace = paste0("http://www.opengis.net/def/crs/EPSG/0/",ref_syst))
  rs$setReferenceSystemIdentifier(rsId)
  md$setReferenceSystemInfo(rs)
  
  # Data identification ----------
  ident <- ISODataIdentification$new()
  ident$setAbstract(product_abstract)
  ident$setPurpose(product_purpose)
  for (f in 1:length(product_credits)){
    ident$addCredit(product_credits[f])
  }
  
  # ident$addCredit("credit2")
  # ident$addCredit("credit3")
  ident$addStatus("completed")
  ident$setLanguage(product_language_code)
  ident$setCharacterSet(product_character_set)
  for (f in 1:length(product_topic_category_codes)) {
    ident$addTopicCategory(product_topic_category_codes[f])
  }
  
  # Point of contact ---------- 
  rp <- ISOResponsibleParty$new()
  # rp$setIndividualName("Kundsupport")
  rp$setOrganisationName(product_org)
  rp$setPositionName("Project Leader")
  rp$setRole("Gustav Kågestan")
  contact <- ISOContact$new()
  phone <- ISOTelephone$new()
  phone$setVoice(product_phone)
  # phone$setFacsimile("myfacsimile")
  contact$setPhone(phone)
  address <- ISOAddress$new()
  address$setDeliveryPoint(product_postal1)
  address$setCity(product_city)
  address$setPostalCode(product_postal2)
  address$setCountry(product_country)
  address$setEmail(product_email)
  contact$setAddress(address)
  res <- ISOOnlineResource$new()
  res$setLinkage(product_base_url)
  res$setName(product_title)
  contact$setOnlineResource(res)
  rp$setContactInfo(contact)
  ident$addPointOfContact(rp)
  
  # Citation ----------
  
  ct <- ISOCitation$new()
  #ct$setTitle(paste0(friendly_title," (", str_split_fixed(pred_name, "_",2)[1],") ",model_area)) ##use file name as title
  ct$setTitle(pred_name1)
  d <- ISODate$new()
  d$setDate(ISOdate(c_yr, c_month, c_day, 1))
  d$setDateType("publication")
  ct$addDate(d)
  ct$setEdition(product_version)
  # ct$setEditionDate(as.Date(ISOdate(c_yr, c_month, c_day, 1)))
  ct$setEditionDate(ISOdate(c_yr, c_month, c_day, 1))
  # ct$setEditionDate(ISOdate(c_yr, c_month, c_day))
  ct$addIdentifier(ISOMetaIdentifier$new(code = "identifier"))
  ct$addPresentationForm("mapDigital")
  ct$addCitedResponsibleParty(rp)
  ident$setCitation(ct)
  
  # Graphic overview ----------
  
  # go1 <- ISOBrowseGraphic$new(
  #   fileName = "http://wwww.somefile.org/png1",
  #   fileDescription = "Map Overview 1",
  #   fileType = "image/png"
  # )
  
  # go2 <- ISOBrowseGraphic$new(
  #   fileName = "http://www.somefile.org/png2",
  #   fileDescription = "Map Overview 2",
  #   fileType = "image/png"
  # )
  # ident$addGraphicOverview(go1)
  # ident$addGraphicOverview(go2)
  
  # Maintenance information ----------
  mi <- ISOMaintenanceInformation$new()
  mi$setMaintenanceFrequency(product_update_frequency)
  ident$setResourceMaintenance(mi)
  
  # Access legal constraints ----------
  # for INSPIRE controlled terms on access legal constraints, please browse the INSPIRE registry:
  # http://inspire.ec.europa.eu/metadata-codelist/LimitationsOnPublicAccess/
  lc <- ISOLegalConstraints$new()
  lc$addAccessConstraint(product_access_constraints)
  
  # lc$addAccessConstraint("otherRestrictions")
  # lc$addOtherConstraint(ISOAnchor$new(
  #   href = "http://inspire.ec.europa.eu/metadata-codelist/LimitationsOnPublicAccess/INSPIRE_Directive_Article13_1a",
  #   name = "public access limited according to Article 13(1)(a) of the INSPIRE Directive"
  # ))
  ident$addResourceConstraints(lc)
  
  # adding use legal constraints
  # for INSPIRE controlled terms on use legal constraints, please browse the INSPIRE registry:
  # http://inspire.ec.europa.eu/metadata-codelist/ConditionsApplyingToAccessAndUse
  lc2 <- ISOLegalConstraints$new()
  for (f in 1:length(product_use_limitation)) {
    lc2$addUseLimitation(product_use_limitation[f])
  }
  
  lc2$addAccessConstraint(product_access_constraint)
  lc2$addOtherConstraint(ISOAnchor$new(
    href = product_legal_constraints_url,
    name = product_legal_constraints_desc
  ))
  ident$addResourceConstraints(lc2)
  
  # Security constraints ----------
  
  sc <- ISOSecurityConstraints$new()
  sc$setClassification(product_secrecy)
  sc$setUserNote(product_secrecy_notes)
  sc$setClassificationSystem(product_secrecy_classification)
  sc$setHandlingDescription(product_secrecy_handling)
  ident$addResourceConstraints(sc)
  
  # Extent ----------
  
  if (product_data_type == "grid") {
    # -----Convert to WGS84, won´t work
    x <- c(extent(aRaster)[1], extent(aRaster)[2])
    y <- c(extent(aRaster)[3], extent(aRaster)[4])
    d <- data.frame(lon = x, lat = y)
    coordinates(d) <- c("lon", "lat")
    proj4string(d) <- CRS(curr_ref_syst)
    CRS.new <- CRS("+init=epsg:4326")
    d.new <- spTransform(d, CRS.new)
    iso_extent <- ISOExtent$new()
    min_x <- round(xmin(d.new), 4)
    max_x <- round(xmax(d.new), 4)
    min_y <- round(ymin(d.new), 4)
    max_y <- round(ymax(d.new), 4)
    bbox <- ISOGeographicBoundingBox$new(minx = min_x, maxx = max_x, miny = min_y, maxy = max_y)
    # bbox <- ISOGeographicBoundingBox$new(minx = xmin(d.new), miny =ymin(d.new), maxx = max_x, maxy=ymax(d.new))
    
    # ------ Use existing EPSG, won´t work
    # bbox <- ISOGeographicBoundingBox$new(minx = extent(aRaster)[1], maxx = extent(aRaster)[2], miny = extent(aRaster)[3], maxy = extent(aRaster)[4])
    
    # ------ dummy value, works
    # bbox <- ISOGeographicBoundingBox$new(minx = -180, miny = -90, maxx = 180, maxy = 90)
    
    iso_extent$setGeographicElement(bbox)
    ident$setExtent(iso_extent)
  }
  
  # Keywords ----------
  
  kwds <- ISOKeywords$new()
  for (f in 1:length(product_theme_keywords)) {
    kwds$addKeyword(product_theme_keywords[f])
  }
  kwds$setKeywordType("theme")
  th <- ISOCitation$new()
  th$setTitle("General")
  d <- ISODate$new()
  d$setDate(ISOdate(c_yr, c_month, c_day, 1))
  d$setDateType("publication")
  th$addDate(d)
  kwds$setThesaurusName(th)
  ident$addKeywords(kwds)
  
  # supplementalInformation
  ident$setSupplementalInformation(product_supplemental_info)
  
  # spatial representation type
  ident$addSpatialRepresentationType(product_data_type)
  
  md$setIdentificationInfo(ident)
  
  # Distribution ----------
  
  distrib <- ISODistribution$new()
  dto <- ISODigitalTransferOptions$new()
  or <- ISOOnlineResource$new()
  or$setLinkage(product_distribution_link)
  or$setName(paste0(product_distribution_name))
  or$setDescription(product_distribution_description)
  or$setProtocol(product_distribution_protocol)
  dto$addOnlineResource(or)
  
  distrib$setDigitalTransferOptions(dto)
  md$setDistributionInfo(distrib)
  
  # Data quality & conformance ----------
  
  # create dataQuality object with a 'dataset' scope
  dq <- ISODataQuality$new()
  scope <- ISOScope$new()
  scope$setLevel("dataset")
  dq$setScope(scope)
  
  # add a report the data quality
  dc <- ISODomainConsistency$new()
  result <- ISOConformanceResult$new()
  spec <- ISOCitation$new()
  spec$setTitle(product_data_quality_check)
  spec$addAlternateTitle(product_data_quality_check_report) 
  d <- ISODate$new()
  # d$setDate(as.Date(ISOdate(c_yr, c_month, c_day, 1)))
  d$setDate(ISOdate(c_yr, c_month, c_day, 1))
  d$setDateType("publication")
  spec$addDate(d)
  result$setSpecification(spec)
  result$setExplanation(product_quality_conformance)
  result$setPass(TRUE)
  dc$addResult(result)
  dq$addReport(dc)
  
  # INSPIRE ----------
  
  # add INSPIRE reports?
  # INSPIRE - interoperability of spatial data sets and services
  dc_inspire1 <- ISODomainConsistency$new()
  cr_inspire1 <- ISOConformanceResult$new()
  cr_inspire_spec1 <- ISOCitation$new()
  cr_inspire_spec1$setTitle(inspire_spec1_title)
  cr_inspire1$setExplanation(inspire_explanation)
  cr_inspire_date1 <- ISODate$new()
  cr_inspire_date1$setDate(as.Date(ISOdate(2010, 12, 8)))
  cr_inspire_date1$setDateType("publication")
  cr_inspire_spec1$addDate(cr_inspire_date1)
  cr_inspire1$setSpecification(cr_inspire_spec1)
  cr_inspire1$setPass(TRUE)
  dc_inspire1$addResult(cr_inspire1)
  dq$addReport(dc_inspire1)
  # INSPIRE - metadata
  dc_inspire2 <- ISODomainConsistency$new()
  cr_inspire2 <- ISOConformanceResult$new()
  cr_inspire_spec2 <- ISOCitation$new()
  cr_inspire_spec2$setTitle(inspire_spec2_title)
  cr_inspire2$setExplanation(inspire_explanation)
  cr_inspire_date2 <- ISODate$new()
  cr_inspire_date2$setDate(as.Date(ISOdate(2008, 12, 4)))
  cr_inspire_date2$setDateType("publication")
  cr_inspire_spec2$addDate(cr_inspire_date2)
  cr_inspire2$setSpecification(cr_inspire_spec2)
  cr_inspire2$setPass(TRUE)
  dc_inspire2$addResult(cr_inspire2)
  dq$addReport(dc_inspire2)
  
  # Lineage ----------
  # add lineage (more example of lineages in ISOLineage documentation)
  lineage <- ISOLineage$new()
  lineage$setStatement(product_lineage_statement)
  dq$setLineage(lineage)
  
  md$setDataQualityInfo(dq)
  
  # XML representation of the ISOMetadata
  
  # --- Testa att kommentera bort!!
  #xml <- md$encode()
  
  
  
  md$save(xml_file_path)#, inspire = TRUE)
  
  
  #xml_file_path = 
  #md$save(x)#, inspire = TRUE)
  #xml_file_path <- paste0(x, ".xml"
  # md$encode(inspire = TRUE)
  # md$save("mymetadata.xml", inspire = TRUE)
}


