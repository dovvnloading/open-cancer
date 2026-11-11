"use strict";
// Compact, lossless encoding of public metadata health observations for GitHub
// Issues (65,536-character body limit). Not for sensitive or patient-level data.
const zlib = require("node:zlib");
const MAX_RAW_BYTES = 4 * 1024 * 1024;
const MARKER = "catalog-health-gzip-base64";

function encodeLedger(report) {
  const json = Buffer.from(JSON.stringify(report), "utf8");
  if (json.length > MAX_RAW_BYTES) throw new Error("Catalog health ledger exceeds raw size limit");
  const packed = zlib.gzipSync(json, {level: 9}).toString("base64");
  return "<!-- " + MARKER + "\n" + packed + "\n-->";
}

function decodeLedger(body) {
  if (typeof body !== "string") return null;
  const compressed = body.match(/<!-- catalog-health-gzip-base64\n([A-Za-z0-9+/=]+)\n-->/);
  if (compressed) {
    const encoded = compressed[1];
    if (encoded.length > MAX_RAW_BYTES * 2) throw new Error("Compressed ledger oversized");
    const content = zlib.gunzipSync(Buffer.from(encoded, "base64"), {
      maxOutputLength: MAX_RAW_BYTES,
    });
    return JSON.parse(content.toString("utf8"));
  }
  // Backward compatibility with older ledger issues created before compression.
  const legacy = body.match(/<!-- catalog-health-json\n([\s\S]*?)\n-->/);
  return legacy ? JSON.parse(legacy[1]) : null;
}

module.exports = {encodeLedger, decodeLedger};
