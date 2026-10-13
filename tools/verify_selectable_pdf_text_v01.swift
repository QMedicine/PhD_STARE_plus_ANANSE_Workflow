// Laptop document validation only; scientific plots are generated on the HPC.
import Foundation
import PDFKit
let root=URL(fileURLWithPath:CommandLine.arguments[1])
let receiptData=try Data(contentsOf:root.appendingPathComponent("selected_report_qualification.json"))
let receipt=try JSONSerialization.jsonObject(with:receiptData) as! [String:Any]
let expected=receipt["pdf_count"] as! Int
var warningFiles=[String:(String,Double)]()
for record in receipt["records"] as! [[String:Any]] {
 let parent=((record["report"] as! String as NSString).deletingLastPathComponent as NSString).deletingLastPathComponent
 let settings=record["settings"] as! [String:Any]
 for tf in record["center_only_TFs"] as! [String] {
  let name="ANANSE_2step_branch_\(tf)_\(record["source"]!)_to_\(record["target"]!).pdf"
  warningFiles[parent+"/plots/"+name]=(tf,(settings["min_edge_weight_plot"] as! NSNumber).doubleValue)
 }
}
var records=[[String:Any]]()
let files=FileManager.default.enumerator(at:root,includingPropertiesForKeys:nil)!
for case let url as URL in files where url.pathExtension=="pdf" {
 guard let doc=PDFDocument(url:url),let contents=doc.string,doc.pageCount>0,contents.count>50 else { fatalError("Invalid or non-selectable PDF: \(url.path)") }
 let normalized=contents.precomposedStringWithCompatibilityMapping
 let text=normalized.components(separatedBy:.whitespacesAndNewlines).filter{!$0.isEmpty}.joined(separator:" ")
 let relative=url.path.replacingOccurrences(of:root.path+"/",with:"")
 if let (tf,cutoff)=warningFiles[relative] {
  let phrase="Warning, \(tf) edge values fall below the minimum cutoff value of "
  guard let range=text.range(of:phrase),let actual=Double(text[range.upperBound...].split(separator:" ")[0]),abs(actual-cutoff)<1e-12 else { fatalError("Warning text or threshold missing: \(url.path)") }
 }
 records.append(["file":relative,"pages":doc.pageCount,"selectable_text_characters":contents.count,"warning":warningFiles[relative] != nil])
}
guard records.count==expected else { fatalError("Expected \(expected) PDFs, found \(records.count)") }
let data=try JSONSerialization.data(withJSONObject:["status":"PASS_SELECTABLE_PDF_TEXT","pdf_count":records.count,"verified_warning_panels":warningFiles.count,"records":records],options:[.prettyPrinted,.sortedKeys])
try data.write(to:root.appendingPathComponent("pdf_text_validation.json"))
print("Verified \(records.count) selectable PDFs and \(warningFiles.count) exact warning panels")
